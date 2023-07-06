
#ifndef ATTENUATION_TESTS_H
#define ATTENUATION_TESTS_H

#include "gtest/gtest.h"

#include "BinaryFileReader.h"

#include <memory>

class AttenuationTests : public testing::Test {
protected:
	AttenuationTests() {}

	// Per-test-suite set-up.
	// Called before the first test in this test suite.
	// Can be omitted if not needed.
	static void SetUpTestCase() {
		// Avoid reallocating static objects if called in subclasses of AttenuationTests.
		if (m_ituDataFileReader == nullptr) {
			// NOTE: Update this directory to the location of the itudata\ directory containing all of the BINGRD data grid files required for ITU model computations
			const std::string filePath = "/home/ayeh/itu/ituModels/itu_excel_example/itudata";
			//"C:\\Users\\Joel Ortiz\\source\\repos\\ituModels\\itu_excel_example\\itudata";
			m_ituDataFileReader = std::make_shared<BinaryFileReader>(filePath);
			std::cout << "INFO: AttenuationTests::AttenuationTests(): Finished loading the BinaryFileReader pointer to the itudata/ directory!" << std::endl;
		}
	}

	// Per-test-suite tear-down.
	// Called after the last test in this test suite.
	// Can be omitted if not needed.
	static void TearDownTestCase() { }

	// You can define per-test set-up logic as usual.
	void SetUp() override { }

	// You can define per-test tear-down logic as usual.
	void TearDown() override { }

	// Data management object that should only be initialized once for all unit tests
	static std::shared_ptr<BinaryFileReader> m_ituDataFileReader;
};

#endif /* ATTENUATION_TESTS_H */