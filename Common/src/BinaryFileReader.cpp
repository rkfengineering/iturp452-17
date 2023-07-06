#include "BinaryFileReader.h"

#include <bit>
#include <cstdint>
#include <climits>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

BinaryFileReader::BinaryFileReader(const std::string & parentDirPath) : _parentDirectory(parentDirPath), _isLittleEndian(false) {
    // Big endian means that the most significant bit comes first
    if constexpr (std::endian::native == std::endian::big) {
        _isLittleEndian = false;
    }
    // Little endian means that the least significant bit comes first
    else if constexpr (std::endian::native == std::endian::little) {
        _isLittleEndian = true;
    }
}

template<typename T>
T BinaryFileReader::readBinaryData(std::istream &inStream, const bool &needSwapForEndianness) const {
    // Make sure that a character takes up 8 bits (1 byte) on this system
    static_assert (CHAR_BIT == 8, "CHAR_BIT != 8");

    const size_t typeSize_bytes = sizeof(T);

    // Use a union to aggregate individual 1 byte values into one larger final value
    union{
        T finalValue; // This will be size of T (which has a size of N bytes)
        uint8_t oneByteArray[typeSize_bytes]; // This will have an array of length N with 1 byte values
    } unionAggregator;

    // Extract contents of binary file up until we've extracted # bits = size of T
    for (uint32_t oneByteInd = 0; oneByteInd < typeSize_bytes; oneByteInd++){
        unionAggregator.oneByteArray[oneByteInd] = inStream.get();
    }

    // If necessary, swap the bit order of the oneByteArray if endianness on this system is backwards (i.e., little endian)
    if(needSwapForEndianness){
        for (uint32_t oneByteInd = 0; oneByteInd < (typeSize_bytes / 2); oneByteInd++) {
            // Swap first value in the array with the last value in the array, 
            // clamping down towards the middle until the first half of values have been swapped with the latter half of values
            std::swap(
                unionAggregator.oneByteArray[oneByteInd],
                unionAggregator.oneByteArray[typeSize_bytes - 1 - oneByteInd]
            );
        }
    }

    // Now that we've aggregated together a series of 1-byte values, we can extract our final N-byte value
    return unionAggregator.finalValue;
}

std::vector<std::vector<double>> BinaryFileReader::readGridData(const std::string& fileName) const {
    // Concatenate directory and filename into full file path
    std::filesystem::path directory(_parentDirectory);
    std::filesystem::path file(fileName + ".bingrd");
    std::filesystem::path filePath = directory / file;

    std::ifstream fileStream;
    // Attempt to open the binary file
    fileStream.exceptions(fileStream.exceptions() | std::ios::failbit);
    try {
        fileStream.open(filePath, std::ios::in | std::ios::binary);
    }
    catch (std::exception& err) {
        std::ostringstream oStrStream;
        oStrStream << "ERROR: BinaryFileReader::readGridData(): Failed to open file in \"" << filePath << "\": " << err.what() << "!";
        throw std::runtime_error(oStrStream.str());
    }
    // Determine size of the file
    fileStream.seekg(0, std::ios::end);
    const std::streampos fileSize = fileStream.tellg();
    // Reset to the beginning of the file
    fileStream.seekg(0, std::ios::beg);
    // Assuming the first 32 bits of the file contains the number of rows, get the row count
    uint32_t rowCount = readBinaryData<uint32_t>(fileStream, _isLittleEndian);
    // Assuming the second 32 bits of the file contains the number of columns, get the column count
    uint32_t colCount = readBinaryData<uint32_t>(fileStream, _isLittleEndian);

    double defaultValue = 0.0;

    std::vector<std::vector<double>> fileContents(rowCount);
    // Now, process the remaining data in the file as the "contents"
    for (uint32_t rowInd = 0; rowInd < rowCount; ++rowInd) {
        // Initialize current row's vector to a vector of size colCount with default values of 0
        fileContents[rowInd].resize(colCount, defaultValue);
        for (uint32_t colInd = 0; colInd < colCount; ++colInd) {
            fileContents[rowInd][colInd] = readBinaryData<double>(fileStream, _isLittleEndian);
        }
    }

    // If we haven't reached the end of this binary file, something probably went wrong
    if (fileStream.tellg() != fileSize) {
        throw std::runtime_error("ERROR: BinaryFileReader::readGridData(): Detected more data in the binary file than expected!");
    }

    return fileContents;
}
