#pragma once

#include <ryml.hpp>
#include <ryml_std.hpp>

namespace OTAP
{
    template <class CharContainer = std::string>
    size_t file_get_contents(const char *filename, CharContainer *v)
    {
        std::FILE *fp = std::fopen(filename, "rb");
        RYML_CHECK_MSG(fp != nullptr, "could not open file");
        std::fseek(fp, 0, SEEK_END);
        long sz = std::ftell(fp);
        v->resize(static_cast<typename CharContainer::size_type>(sz));
        if (sz)
        {
            std::rewind(fp);
            size_t ret = std::fread(&(*v)[0], 1, v->size(), fp);
            RYML_CHECK(ret == (size_t)sz);
        }
        std::fclose(fp);
        return v->size();
    }

    /** load a file from disk and return a newly created CharContainer */
    template <class CharContainer = std::string>
    CharContainer file_get_contents(const char *filename)
    {
        CharContainer cc;
        file_get_contents(filename, &cc);
        return cc;
    }

    /** save a buffer into a file */
    template <class CharContainer = std::string>
    static void file_put_contents(const char *filename, CharContainer const &v, const char *access)
    {
        file_put_contents(filename, v.empty() ? "" : &v[0], v.size(), access);
    }

    /** save a buffer into a file */
    void file_put_contents(const char *filename, const char *buf, size_t sz, const char *access)
    {
        std::FILE *fp = std::fopen(filename, access);
        RYML_CHECK_MSG(fp != nullptr, "could not open file");
        std::fwrite(buf, 1, sz, fp);
        std::fclose(fp);
    }
}