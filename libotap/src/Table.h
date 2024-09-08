#pragma once

#include <vector>
#include "rapidcsv/rapidcsv.h"
#include "mlinterp/mlinterp.hpp"

namespace OTAP
{
    // NOTE: May extend to multiple columns with different datatypes
    // TEST: Check edge cases and test this
    template <typename T, size_t dim = 2>
    class Table
    {
        // using Type = std::enable_if_t<std::is_arithmetic_v<T>, T>;
        using Type = std::remove_reference_t<std::remove_cv_t<T>>;
        static_assert(std::is_arithmetic_v<T>, "Table must be of arithmetic type");
        static_assert(dim >= 2, "Number of columns should be atleast two");
        using columnType = std::vector<Type>;

        std::vector<columnType> mData;
        std::map<Type, size_t> mColumnNames;
        std::map<Type, size_t> mRowNames;

    public:
        // Table() = default;
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

        ~Table() = default;
    };
} // namespace OTAP
