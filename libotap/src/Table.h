#pragma once

#include <vector>
#include "Eigen/Eigen"
#include "rapidcsv/rapidcsv.h"
#include "mlinterp/mlinterp.hpp"

namespace OTAP
{
    // NOTE: May extend to multiple columns with different datatypes
    // TEST: Check edge cases and test this
    // FIXME: Add support for rutime dynamic sizing
    template <typename T = double, size_t dim = 2, typename C = size_t>
    class Table
    {
        // using Type = std::enable_if_t<std::is_arithmetic_v<T>, T>;
        using Type = std::remove_reference_t<std::remove_cv_t<T>>;
        static_assert(std::is_arithmetic_v<T>, "Table must be of arithmetic type");
        static_assert(dim >= 2, "Number of columns should be atleast two");
        using columnType = std::vector<Type>;

        std::vector<columnType> mData;
        std::map<Type, C> mColumnNames;
        std::map<Type, C> mRowNames;

    public:
        explicit Table() = default;
        Table(const std::string &filepath, bool delim_whitespace = true)
        {
            Load(filepath, delim_whitespace);
        }

        // TODO: Chek validity with brace-enclosed initializer lists
        Table(std::initializer_list<columnType> args)
        {
            mData = args;
        }

        Table(std::istream &stream, bool delim_whitespace = true)
        {
            Load(stream, delim_whitespace);
        }

        void Load(const std::string &filepath, bool delim_whitespace = true)
        {
            rapidcsv::Document doc(filepath, rapidcsv::LabelParams(-1, -1),
                                   delim_whitespace ? rapidcsv::SeparatorParams() : rapidcsv::SeparatorParams(false, ',', true));
            for (size_t i = 0; i < doc.GetColumnCount(); i++)
            {
                // mData.push_back(std::stod(doc.GetColumnName(i))); // FIXME: always passes double
                mData.push_back(doc.GetColumn<Type>(i));
            }
        }
        void Load(std::istream &stream, bool delim_whitespace = true)
        {
            rapidcsv::Document doc(stream, rapidcsv::LabelParams(-1, -1),
                                   delim_whitespace ? rapidcsv::SeparatorParams() : rapidcsv::SeparatorParams(false, ',', true));
            for (size_t i = 0; i < doc.GetColumnCount(); i++)
            {
                // mData.push_back(std::stod(doc.GetColumnName(i))); // FIXME: always passes double
                mData.push_back(doc.GetColumn<Type>(i));
            }
        }
        // FIXME: Add saving functions for Table
        // void Save(const std::string &filepath, bool delim_whitespace = true)
        // {
        //     rapidcsv::Document doc(filepath, rapidcsv::LabelParams(),
        //                            delim_whitespace ? rapidcsv::SeparatorParams() : rapidcsv::SeparatorParams(false, ',', true));
        //     doc.Save(filepath);
        // }
        // void Save(std::istream &stream, bool delim_whitespace = true)
        // {
        //     rapidcsv::Document doc(stream, rapidcsv::LabelParams(),
        //                            delim_whitespace ? rapidcsv::SeparatorParams() : rapidcsv::SeparatorParams(false, ',', true));
        //     doc.Save(stream);
        // }

        std::conditional_t<dim == 2, Type, std::vector<Type>> operator[](const Type &index)
        {
            if constexpr (dim == 2)
            {
                Type val;
                size_t dims[] = {mData[0].size()};
                mlinterp::interp(dims, size_t(1), mData[1].data(), &val, mData[0].data(), &index);
                return val;
            }
            else
            {
                Type val;
                std::vector<Type> vals;
                size_t dims[] = {mData[0].size()};
                for (size_t i = 1; i < dim; i++)
                {
                    mlinterp::interp(dims, size_t(1), mData[i].data(), &val, mData[0].data(), &index);
                    vals.push_back(val);
                }
                return vals;
            }
        }

        Table<T, dim> &operator=(const std::initializer_list<columnType> &args)
        {
            mData = args;
            return *this;
        }

        ~Table() = default;
    };

    // Table Grid for runtime dynamic tables

    enum class TableGridOrder
    {
        ColumnFirst,
        RowFirst
    };

    template <typename T = double, TableGridOrder O = TableGridOrder::ColumnFirst, TableGridOrder I = TableGridOrder::RowFirst>
    class TableGrid
    {
        using Type = std::remove_reference_t<std::remove_cv_t<T>>;
        using IndexType = size_t;
        using columnType = std::vector<Type>;
        static_assert(std::is_arithmetic_v<T>, "TableGrid must be of arithmetic type");

        std::vector<columnType> mData;
        std::map<Type, IndexType> mColumnNames;
        std::map<Type, IndexType> mRowNames;

        size_t m_nrows = 0;
        size_t m_ncols = 0;

    public:
        TableGrid() = default;
        TableGrid(const std::vector<columnType> &args) : mData(args) {}
        TableGrid(const std::string &filepath, bool delim_whitespace = true) { Load(filepath, delim_whitespace); }
        TableGrid(std::istream &stream, bool delim_whitespace = true) { Load(stream, delim_whitespace); }

        void Load(const std::string &filepath, bool delim_whitespace = true)
        {
            std::ifstream inputfile(filepath);
            Load(inputfile, delim_whitespace);
        }

        void Load(std::istream &stream, bool delim_whitespace = true)
        {
            rapidcsv::Document doc(stream, rapidcsv::LabelParams(-1, -1),
                                   delim_whitespace ? rapidcsv::SeparatorParams() : rapidcsv::SeparatorParams(false, ',', true));
            for (size_t i = 0; i < doc.GetColumnCount(); i++)
            {
                mData.push_back(doc.GetColumn<Type>(i));
            }
            m_nrows = doc.GetRowCount();
            m_ncols = doc.GetColumnCount();
        }

        void Save(const std::string &filepath, bool delim_whitespace = true) const
        {
            std::ofstream outfile(filepath);
            Save(outfile, delim_whitespace);
        }

        void Save(std::ostream &stream, bool delim_whitespace = true) const
        {
            Eigen::MatrixXd mat;
            mat.resize(m_nrows, m_ncols);
            Eigen::IOFormat fmt;
            fmt.matSuffix = "\n";
            if (delim_whitespace)
                fmt.coeffSeparator = "\t";
            else
                fmt.coeffSeparator = ",\t";
            for (size_t i = 0; i < m_ncols; i++)
                for (size_t j = 0; j < m_nrows; j++)
                    mat(j, i) = mData[i][j];
            if constexpr (O == TableGridOrder::ColumnFirst)
                stream << mat.format(fmt);
            else
                stream << mat.transpose().format(fmt);
        }

        columnType operator[](const Type &index) const
        {
            Type val;
            columnType vals;
            if constexpr (I == TableGridOrder::RowFirst)
            {
                size_t dims[] = {m_nrows};
                for (size_t i = 1; i < m_ncols; i++)
                {
                    mlinterp::interp(dims, size_t(1), mData[i].data(), &val, mData[0].data(), &index);
                    vals.push_back(val);
                }
                return vals;
            }
            else
            {
                static_assert(false);
                size_t dims[] = {m_nrows, m_ncols};
                for (size_t i = 1; i < m_nrows; i++)
                {
                    mlinterp::interp(dims, m_nrows, mData.data(), &val, mData[0].data(), &index);
                    vals.push_back(val);
                }
                return vals;
            }
        }

        Type &operator()(size_t ci, size_t ri) { return mData.at(ci).at(ri); }

        TableGrid<T, O, I> &operator=(const std::vector<columnType> &args)
        {
            mData = args;
            m_ncols = mData.size();
            m_nrows = mData[0].size();
            return *this;
        }

        TableGrid<T, O, I> &addColumn(columnType c)
        {
            assert(c.size() == m_nrows || m_nrows == 0);
            mData.emplace_back(c);
            m_ncols++;
            m_nrows = c.size();
            return *this;
        }

        TableGrid<T, O, I> &insertColumn(IndexType i, columnType c)
        {
            assert(c.size() == m_nrows || m_nrows == 0);
            mData.insert(mData.begin() + i, c);
            m_ncols++;
            m_nrows = c.size();
            return *this;
        }

        TableGrid<T, O, I> &deleteColumn(IndexType i)
        {
            mData.erase(mData.begin() + i);
            m_ncols--;
            assert(m_ncols >= 2);
            return *this;
        }

        TableGrid<T, O, I> &addRow(columnType c)
        {
            assert(c.size() == m_ncols || m_ncols == 0);
            for (size_t i = 0; i < m_ncols; i++)
                mData[i].emplace_back(c[i]);
            m_nrows++;
            m_ncols = c.size();
            return *this;
        }

        TableGrid<T, O, I> &insertRow(IndexType i, columnType c)
        {
            assert(c.size() == m_ncols || m_ncols == 0);
            for (size_t j = 0; j < m_ncols; j++)
                mData[j].insert(mData[j].begin() + i, c[j]);
            m_nrows++;
            m_ncols = c.size();
            return *this;
        }

        TableGrid<T, O, I> &deleteRow(IndexType i)
        {
            for (size_t j = 0; j < m_ncols; j++)
                mData[j].erase(mData[j].begin() + i);
            m_nrows--;
            assert(m_nrows >= 0);
            return *this;
        }

        const std::vector<columnType> &data() const { return mData; }

        size_t getRowCount() const { return m_nrows; }
        size_t getColumnCount() const { return m_ncols; }

        friend std::ostream &operator<<(std::ostream &stream, const TableGrid<T, O, I> &tab)
        {
            tab.Save(stream);
            return stream;
        }

        ~TableGrid() = default;
    };

} // namespace OTAP
