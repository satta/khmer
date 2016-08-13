// g++ --std=c++11 -O3 -I ../lib -I ../third-party/seqan/core/include/ -o snp-power snp-power.cpp ../lib/liboxli.a

#include <assert.h>
#include <getopt.h>
#include <cmath>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>

#include "hashbits.hh"
#include "khmer.hh"
#include "khmer_exception.hh"
#include "read_parsers.hh"

using namespace khmer;
using namespace read_parsers;

typedef std::vector<HashIntoType> vint;
typedef std::vector<std::string> vstring;
typedef Read Sequence;

class Mutator
{
public:
    Sequence& seq;
    int k;
    int i;
    std::string nucl;

    Mutator(Sequence& sequence, int ksize)
        : seq(sequence), k(ksize), i(k - 1), nucl("ACGT") {}

    bool next_mutated_seqs(vstring& mutseqs)
    {
        mutseqs.clear();
        if (i + k > seq.sequence.length()) {
            return false;
        }

        int minpos = i - k + 1;
        int len = (2*k) - 1;
        std::string window = seq.sequence.substr(minpos, len);

        for (auto bp : nucl) {
            if (window[k - 1] != bp) {
                std::string mutseq = window;
                mutseq[k - 1] = bp;
                mutseqs.push_back(mutseq);
            }
        }

        i++;
        return true;
    }
};

bool isprime(int n)
{
    assert(n >= 0);
    if (n == 1) {
        return false;
    }
    else if (n == 2) {
        return true;
    }
    else if (n % 2 == 0) {
        return true;
    }
    for (int i = 3; i < pow(n, 0.5) + 1; i += 2) {
        if (n % i == 0) {
            return false;
        }
    }
    return true;
}

vint get_n_primes_near_x(int x, int n)
{
    vint primes;
    if (x == 1 && n == 1) {
        primes.push_back(1);
        return primes;
    }

    int i = x - 1;
    if (i % 2 == 0) {
        i--;
    }
    while (primes.size() < n && i > 0) {
        if (isprime(i)) {
            primes.push_back(i);
        }
        i -= 2;
    }
    assert(primes.size() == n);
    return primes;
}

void count_mutation_collisions(Hashbits& nodegraph, IParser *parser, int k,
                               unsigned long limit)
{
    unsigned long long total = 0;
    unsigned long hits = 0;
    Sequence seq;
    unsigned long count = 0;

    while (!parser->is_complete()) {
        try {
            seq = parser->get_next_read();
        } catch (NoMoreReadsAvailable &e) {
            break;
        }

        Mutator m(seq, k);
        vstring mutseqs;
        while (m.next_mutated_seqs(mutseqs) && (limit == 0 || count < limit)) {
            count++;
            if (count % 100000 == 0) {
                std::cerr << "  processed " << (float)count / (float)1000
                          << " Kb of sequence" << std::endl;
            }
            vstring::iterator it;
            for (it = mutseqs.begin(); it < mutseqs.end(); it++) {
                std::string mutseq = *it;
                vstring kmers;
                nodegraph.get_kmers(mutseq, kmers);
                vstring::iterator kit;
                for (kit = kmers.begin(); kit < kmers.end(); kit++) {
                    std::string kmer = *kit;
                    total++;
                    if (nodegraph.get_count(kmer.c_str()) > 0) {
                        hits++;
                    }
                }
            }
        }
    }
    std::cout << total << " " << hits << " " << (float)hits / (float)total
              << std::endl;
}

void print_usage(std::ostream& stream = std::cerr)
{
    stream << "Usage: snp-power [options] seqs.fa [genome.fa]" << std::endl;
    stream << "  options:" << std::endl;
    stream << "    -h    print this help message and exit" << std::endl;
    stream << "    -k    k-mer length (default: 31)" << std::endl;
    stream << "    -x    approx table size (default: 5e8)" << std::endl;
    stream << "    -N    num tables (default: 4)" << std::endl;
    stream << "    -l    limit number of positions to process" << std::endl;
}

int main(int argc, const char **argv)
{
    assert(isprime(37));
    assert(!isprime(51));
    assert(isprime(15485863));
    assert(!isprime(15485865));

    if (argc == 1) {
        print_usage(std::cout);
        return 0;
    }

    int k = 31;
    int targetsize = 5e8;
    int numtables = 4;
    unsigned long limit = 0;

    char c;
    while ((c = getopt (argc, (char *const *)argv, "hk:x:N:l:")) != -1) {
        if (c == 'h') {
            print_usage(std::cout);
            return 0;
        }
        else if (c == 'k') {
            k = atoi(optarg);
        }
        else if (c == 'x') {
            targetsize = atoi(optarg);
        }
        else if (c == 'N') {
            numtables = atoi(optarg);
        }
        else if (c == 'l') {
            limit = atoi(optarg);
        }
        else {
            std::cerr << "Unknown option '" << c << "'" << std::endl;
            print_usage();
            return 1;
        }
    }

    std::string infile(argv[optind]);
    std::string refrfile(argv[optind]);
    if (argc > optind + 1) {
        refrfile = argv[optind + 1];
    }

    std::cerr << "allocating nodegraph" << std::endl;
    vint tablesizes = get_n_primes_near_x(targetsize, numtables);
    Hashbits nodegraph(k, tablesizes);

    std::cerr << "consuming input" << std::endl;
    unsigned int seqs_consumed = 0;
    unsigned long long kmers_consumed = 0;
    nodegraph.consume_fasta(refrfile, seqs_consumed, kmers_consumed);
    std::cerr << "consumed " << seqs_consumed << " sequence(s) and "
              << kmers_consumed << " " << k << "-mers" << std::endl;

    std::cerr << "generating mutations" << std::endl;
    IParser *parser = IParser::get_parser(infile);
    count_mutation_collisions(nodegraph, parser, k, limit);
    delete parser;

    return 0;
}
