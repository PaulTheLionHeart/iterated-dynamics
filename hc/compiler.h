#pragma once

#include "modes.h"

#include <string>
#include <vector>

namespace hc
{

struct Options
{
    modes mode{modes::NONE};
    std::string fname1;
    std::string fname2;
    std::string swappath;
    bool show_mem{};
    bool show_stats{};
    bool quiet_mode{}; // true if "/Q" option used
    std::vector<std::string> include_paths;
    std::string output_dir;
};

Options parse_options(int argc, char **argv);

class compiler
{
public:
    compiler(const Options &options);
    ~compiler();

    int process();

private:
    void read_source_file();
    void usage();
    void compile();
    void print();
    void render_html();
    void paginate_html_document();
    void print_html_document(std::string const &output_filename);

    Options m_options;
};

} // namespace hc
