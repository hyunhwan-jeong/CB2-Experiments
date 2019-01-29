#include <unordered_map>
#include <map>
#include <fstream>
#include <cstdio>
#include <iostream>
#include <string> 
#include <cctype> 
#include <cstring>
#include <algorithm>
#include <cassert>
#include <sys/stat.h>
#include <unistd.h>
#include <libgen.h>
#include <ctime>
#include "fmt/format.h"
#include "fmt/time.h"
using namespace std;

double GetTickCount(void) {
    struct timespec now;
    if (clock_gettime(CLOCK_MONOTONIC, &now))
        return 0;
    return now.tv_sec * 1000.0 + now.tv_nsec / 1000000.0;
}


//https://stackoverflow.com/questions/2336242/recursive-mkdir-system-call-on-unix
int mkpath(const string &str, mode_t mode) {
    char *file_path = new char[str.size()+1];
    strcpy(file_path, str.c_str());
    assert(file_path && *file_path);
    char* p;
    for (p=strchr(file_path+1, '/'); p; p=strchr(p+1, '/')) {
        *p='\0';
        if (mkdir(file_path, mode)==-1) {
            if (errno!=EEXIST) { *p='/'; return -1; }
        }
        *p='/';
    }
    delete [] file_path;
    return 0;
}

inline bool is_exists (const std::string& name) {
    return ( access( name.c_str(), F_OK ) != -1 );
}

string str_dirname(const std::string &file) {
    char* cstr = new char[file.size()+1];
    strcpy(cstr, file.c_str());
    char *out_dir = dirname(cstr);
    string ret = out_dir;
    delete [] cstr;
    return ret;
}

class Logger {
    private:
        ofstream file;
        ostream &oup;
    public:
        Logger(const char *f_log) : file(f_log), oup(file) {}
        Logger() : oup(cerr) {}
        void add_log(const string str) {
            std::time_t t = std::time(nullptr);
            oup << fmt::format("[{:%F %T}] {}", *std::localtime(&t), str) << endl;
        }
};

class sgRNA_MAP {
    private:
        unordered_map<long long, string> lib;
        unordered_map<long long, int> cnt;
        unordered_map<long long, int> cnt_rc;
        unordered_map<int, int> pos;
        unordered_map<int, int> pos_rc;
        int lib_seq_len;
        int num_proc_line;
        int num_hits;
        int num_hits_rc;
        long long tot_reads_len;
        double tot_run_time;
        bool is_rc;

        void load_library(const char *f_lib) {
            ifstream inp(f_lib);
            string name, se;
            lib_seq_len = 20;
            while(inp>>name) {
                if(name[0]=='>') name = name.substr(1);
                inp >> se;
                long long num = 0;
                for(auto &x : se) {
                    num *= 4;
                    num += (toupper(x)>>1)&3;
                }
                lib_seq_len = se.size();
                lib[num] = name;
                cnt[num] = 0;
                cnt_rc[num] = 0;
            }
            inp.close();
        }


        void run_MAP(const char *f_seq, const char *f_log, const char *f_unmap) {
            string line;
            int num_line = 0;
            is_rc = false;
            num_proc_line = 0;
            num_hits = 0;
            num_hits_rc = 0;
            tot_reads_len = 0;
            long long mod = 1LL<<(2*lib_seq_len);
            const int rc[] = {2, 3, 0, 1};
            Logger logger_file(f_log);
            Logger logger_cerr;
            string msg;	
            msg = fmt::format("Detects the length of guided RNA is {}.", lib_seq_len);
            logger_file.add_log(msg);
            logger_cerr.add_log(msg); 

            msg = fmt::format("{} guide RNAs were found.", cnt.size());
            logger_file.add_log(msg);
            logger_cerr.add_log(msg); 

            ifstream inp(f_seq);

            ofstream oup_unmap;
            if(strlen(f_unmap)>0) {
                oup_unmap.open(f_unmap, ofstream::out);
            }
            while(getline(inp, line)) {
                if(num_line++%4!=1) continue;
                tot_reads_len += line.size();
                if(++num_proc_line%int(1e6)==0) {
                    msg = fmt::format("Processing {} lines...", num_proc_line);
                    logger_file.add_log(msg);
                    logger_cerr.add_log(msg); 

                    msg = fmt::format("Current mappability : {:.2f}%", 100.0*max(num_hits,num_hits_rc)/(num_proc_line-1));
                    logger_file.add_log(msg);
                    logger_cerr.add_log(msg); 
                }
                int cur_len = 0;
                long long num = 0;
                int i = 0;
                bool is_found = false;
                for(auto &c: line) {
                    i++; 
                    c = toupper(c);
                    if(c=='N') {
                        cur_len = 0;
                        num = 0;
                        continue;
                    }
                    num *= 4;
                    num += (c>>1)&3;
                    num %= mod;
                    if(++cur_len==lib_seq_len) {
                        if(cnt.count(num)>0) {
                            pos[i-lib_seq_len]++; 
                            ++num_hits;
                            cnt[num]++;
                            is_found = true;  
                            break;
                        }
                        --cur_len;
                    }
                }
                is_found = false;
                cur_len = 0;
                num = 0;
                reverse(line.begin(), line.end());
                i = 0;
                for(auto &c: line) {
                    c = toupper(c);
                    if(c=='N') {
                        cur_len = 0;
                        num = 0;
                        continue;
                    }
                    num *= 4;
                    num += rc[(int(c)>>1)&3];
                    num %= mod;
                    if(++cur_len==lib_seq_len) {
                        if(cnt_rc.count(num)>0) {
                            pos_rc[i]++;
                            ++num_hits_rc;
                            cnt_rc[num]++;
                            is_found = true;
                            break;
                        }
                        --cur_len;
                    }
                    i++;
                }
                reverse(line.begin(), line.end());
                if(!is_found) oup_unmap << line << endl;
            }
            msg = fmt::format("All {} reads were proceed!", num_proc_line);
            logger_file.add_log(msg);
            logger_cerr.add_log(msg); 
            if(oup_unmap.is_open()) oup_unmap.close();
            inp.close();

            if(num_hits < num_hits_rc) {
                swap(num_hits, num_hits_rc);
                swap(cnt, cnt_rc);
                swap(pos, pos_rc);
                is_rc = true;
            }
        }
    public:
        sgRNA_MAP(const char *f_lib, const char *f_seq, const char *f_log, const char *f_unmap) {
            double st = GetTickCount();
            load_library(f_lib);
            run_MAP(f_seq, f_log, f_unmap);
            double ed = GetTickCount();
            tot_run_time = (ed-st) / 1000.0;
        }

        void print_count(string file_count) {
            ofstream oup(file_count.c_str());
            vector< pair<string, int> > vec;
            vec.reserve(cnt.size());
            for(auto &x: cnt) {
                vec.push_back(make_pair(lib[x.first], x.second));
            }
            sort(vec.begin(), vec.end());

            for(auto &x: vec) {
                oup << fmt::format("{}\t{}", x.first, x.second) << endl;
            }
            oup.close();
        }

        void print_summary(string file_summary) {
            ofstream oup(file_summary.c_str());
            oup << "Total gRNA counts" << "\t" << lib.size() << endl; 
            oup << "gRNA BP" << "\t" << lib_seq_len << endl;
            oup << "average reads BP" << "\t" << tot_reads_len * 1. / num_proc_line << endl;
            oup << "Total reads count" << "\t" << num_proc_line << endl;
            oup << "Mapping success rate" << "\t" << (num_hits*100.0/num_proc_line) << endl;
            oup << "Total mapped reads" << "\t" << num_hits << endl; oup << "gRNA coverage" << "\t" << num_hits * 1. / lib.size() << endl;
            oup << "Running time" << "\t" << tot_run_time << endl;
            oup << "Seconds for millon reads" << "\t" << tot_run_time / num_proc_line * 1e6 << endl;
            oup << "reverse-complement" << "\t" << is_rc << endl;
            oup << "Total alterative reads count" << "\t" << num_hits_rc << endl;
            oup.close(); 
        }

        void print_mpos(string file_pos) {
            ofstream oup(file_pos.c_str());
            vector< pair<int,int> > vec;
            for(auto &x: pos) {
                int p = x.first+1;
                int c = x.second;
                vec.push_back(make_pair(-c, p));
            }
            sort(vec.begin(), vec.end());
            for(auto &x: vec) {
                oup << fmt::format("{}\t{}", x.second, -x.first) << endl;
            }
            oup.close();
        }
};

int main(int argc, char *argv[]) {
    std::ios_base::sync_with_stdio (false);

    map<char, string> opts;
    for( int c ; (c = getopt (argc, argv, "f:l:o:p:u")) != -1 ; ) {
        if(optarg) {
            opts[c] = optarg;
        }
        else {
            opts[c] = "";
        }
    }

    if(opts.count('p')==0) {
        opts['p'] = "quant";
    }

    if('/'!=*opts['o'].rbegin()) opts['s'] += "/";

    opts['c'] = opts['o'] + opts['p'] + "_count.txt";
    opts['s'] = opts['o'] + opts['p'] + "_summary.txt";
    opts['g'] = opts['o'] + opts['p'] + "_log.txt";
    opts['m'] = opts['o'] + opts['p'] + "_mpos.txt";

    if(opts.count('u')>0) {
        opts['u'] = opts['o'] + opts['p'] + "_umapped.txt";
    }
    else {
        opts['u'] = "";
    }
    string out_dir = str_dirname(opts['c']);
    if(!is_exists(out_dir)) {
        mkpath(opts['c'].c_str(), 0755);
    }

    sgRNA_MAP sgrna_map(opts['l'].c_str(), opts['f'].c_str(), opts['g'].c_str(), opts['u'].c_str());
    sgrna_map.print_count(opts['c']);
    sgrna_map.print_summary(opts['s']);
    sgrna_map.print_mpos(opts['m']);
    return 0;
}
