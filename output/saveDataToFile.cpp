#include "saveDataToFile.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <variant>
#include <type_traits>
#include <sstream>
#include <utility>

// datatypes supported by this variant
using ColumnVariant = std::variant<
    const double *,
    const float *,
    const int *,
    const std::string *,
    const bool *>;

struct ValueWriter
{
    std::ostream &os;
    size_t row;
    int precision;
    NumberFormat format;
    size_t width;

    template <typename T>
    void operator()(const T *columnData)
    {
        if constexpr (std::is_floating_point_v<T>)
        {
            switch (format)
            {
            case NumberFormat::FIXED:
                os << std::fixed;
                break;
            case NumberFormat::SCIENTIFIC:
                os << std::scientific;
                break;
            }
            os << std::setprecision(precision);
        }
        os << std::setw(width) << columnData[row];

        if constexpr (std::is_floating_point_v<T>)
        {
            os.unsetf(std::ios_base::floatfield);
        }
    }
};

void saveDataToFile(
    const std::string &filename,
    const ColumnVariant *columnData,
    const std::string *columnNames,
    const int *precisions,
    const NumberFormat *formats,
    const int *columnWidths,
    size_t numColumns,
    size_t numRows,
    bool writeHeader)
{
    // Validate input parameters
    if (filename.empty())
    {
        std::cerr << "Error: Filename cannot be empty\n";
        return;
    }
    if (!columnData || !columnNames || !precisions || !formats || !columnWidths)
    {
        std::cerr << "Error: Null pointer passed as argument\n";
        return;
    }
    if (numColumns == 0 || numRows == 0)
    {
        std::cerr << "Error: Number of columns and rows must be > 0\n";
        return;
    }

    // Open output file
    std::ofstream outFile(filename, writeHeader ? std::ios::out : std::ios::app);
    if (!outFile.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << "\n";
        return;
    }

    // Calculate maximum width for each column
    std::vector<size_t> maxNumWidths(numColumns, 0);
    for (size_t col = 0; col < numColumns; ++col)
    {
        for (size_t row = 0; row < numRows; ++row)
        {
            std::ostringstream oss;
            ValueWriter writer{oss, row, precisions[col], formats[col], 0};
            std::visit(writer, columnData[col]);

            size_t currentWidth = oss.str().length();
            if (currentWidth > maxNumWidths[col])
            {
                maxNumWidths[col] = currentWidth;
            }
        }
    }

    // Write header if requested
    if (writeHeader)
    {
        // Calculate total column widths
        size_t *totalWidths = new size_t[numColumns];
        for (size_t i = 0; i < numColumns; ++i)
        {
            size_t width1 = static_cast<size_t>(columnWidths[i]);
            size_t width2 = maxNumWidths[i];
            size_t width3 = columnNames[i].length();

            totalWidths[i] = std::max(width1, std::max(width2, width3));
        }

        // Write column index header
        outFile << "#  idx   ";

        // Write centered column names
        for (size_t i = 0; i < numColumns; ++i)
        {
            size_t nameLength = columnNames[i].length();
            size_t padding = totalWidths[i] - nameLength;
            size_t leftPad = padding / 2;
            size_t rightPad = padding - leftPad;

            outFile << std::string(leftPad, ' ')
                    << columnNames[i]
                    << std::string(rightPad, ' ')
                    << "    ";
        }
        outFile << "\n";

        // Calculate and write separator line
        size_t totalHeaderWidth = 8; // "#idx    "
        for (size_t i = 0; i < numColumns; ++i)
        {
            totalHeaderWidth += totalWidths[i] + 4; // data + spacing
        }
        outFile << "#" << std::string(totalHeaderWidth, '=') << "\n";
    }

    // Write data rows
    for (size_t row = 0; row < numRows; ++row)
    {
        // Write row index
        outFile << std::setw(5) << row << "    ";

        // Write each column's data
        for (size_t col = 0; col < numColumns; ++col)
        {
            ValueWriter writer{outFile, row, precisions[col], formats[col], maxNumWidths[col]};
            std::visit(writer, columnData[col]);
            outFile << "    ";
        }
        outFile << "\n";
    }

    outFile.close();
}
