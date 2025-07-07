#ifndef SAVE_DATA_TO_FILE_HPP
#define SAVE_DATA_TO_FILE_HPP

#include <variant>
#include <string>
#include <ios> // For std::ios_base::fmtflags
#include <cstddef> // for size_t

/**
 * @brief Saves columnar data to a file with formatted output
 * 
 * @param filename Output file path
 * @param columnData 2D array of data (columns x rows)
 * @param columnNames Array of column names
 * @param precisions Array of floating point precisions
 * @param formats Array of format flags (scientific, fixed, etc.)
 * @param columnWidths Array of minimum column widths
 * @param numColumns Number of columns
 * @param numRows Number of rows
 * @param writeHeader Whether to write header row
 */

enum class NumberFormat {
    FIXED = std::ios_base::fixed,
    SCIENTIFIC = std::ios_base::scientific
    // Add more formats as needed
};

using ColumnVariant = std::variant<
    const double*,
    const float*,
    const int*,
    const std::string*,
    const bool*
>;

void saveDataToFile(
    const std::string& filename,
    const ColumnVariant* columnData,
    const std::string* columnNames,
    const int* precisions,
    const NumberFormat* formats,
    const int* columnWidths,
    size_t numColumns,
    size_t numRows,
    bool writeHeader = true
);

#endif // SAVE_DATA_TO_FILE_HPP