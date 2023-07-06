#ifndef BINARY_FILE_READER_H
#define BINARY_FILE_READER_H

#include <istream>
#include <string>
#include <vector>

/** A default data source which reads files under a single directory.
 * All files must be named "<name>.bingrd" in the directory.
 */
class BinaryFileReader {
public:
    BinaryFileReader(const std::string& parentDirPath);

    std::vector<std::vector<double>> readGridData(const std::string &fileName) const;

private:
    template<typename T> T readBinaryData(std::istream& inStream, const bool& needSwapForEndianness) const;

    /// Directory for file search
    std::string _parentDirectory;
    /// Endianness detection
    bool _isLittleEndian;
};

#endif /* BINARY_FILE_READER_H */
