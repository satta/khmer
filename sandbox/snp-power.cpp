// g++ -Wall --std=c++11 -O3 -I ../lib -I ../third-party/seqan/core/include/ -o snp-power snp-power.cpp ../lib/liboxli.a

#include <assert.h>
#include <getopt.h>
#include <cmath>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>

#include "counting.hh"
#include "khmer.hh"
#include "khmer_exception.hh"
#include "read_parsers.hh"

using namespace khmer;
using namespace read_parsers;

typedef std::vector<int> vint;
typedef std::vector<std::string> vstring;
typedef Read Sequence;

/// A nucleotide that we have mutated in silico, along with the upstream and
/// downstream nucleotides necessary to construct all k-mers overlapping with
/// the mutated nucleotide.
class MutatedKmerInterval
{
public:
    std::string& original;
    std::string mutseq;
    unsigned long unique_kmer_count;
    vint kmer_abund;

    explicit MutatedKmerInterval(const MutatedKmerInterval&) = default;
    MutatedKmerInterval(std::string& sequence, std::string mutatedseq)
        : original(sequence), mutseq(mutatedseq), unique_kmer_count(0) {}

    void check_mutated_kmers(CountingHash& countgraph, bool debug)
    {
        vstring kmers;
        countgraph.get_kmers(mutseq, kmers);
        int nkmers = kmers.size();
        int kmer_offset = 0;
        for (auto kmer : kmers) {
            int kmerfreq = countgraph.get_count(kmer.c_str());
            kmer_abund.push_back(kmerfreq);
            if (kmerfreq == 0) {
                unique_kmer_count++;
            }
            else {
                if (debug) {
                    std::cerr << "DEBUG orig seq: " << original << '\n';
                    std::cerr << "DEBUG mutd seq: " << mutseq << '\n';
                    std::cerr << "DEBUG         : ";
                    for (int i = 0; i < kmer_offset; i++) {
                        std::cerr << ' ';
                    }
                    std::cerr << kmer << "\n\n";
                }
            }
            kmer_offset++;
        }
        assert(kmer_abund.size() == nkmers);
    }
};
typedef std::vector<MutatedKmerInterval> vmut;

/// A nucleotide of interest, along with k - 1 nucleotides on either side; in
/// other words, the interval containing all k-mers overlapping with the
/// nucleotide of interest.
class KmerInterval
{
public:
    std::string& sequence;
    int k;
    std::string nucl;
    vmut mutseqs;

    KmerInterval(std::string& seq, int k) : sequence(seq), k(k), nucl("ACGT")
    {
        for (auto bp : nucl) {
            if (bp != sequence[k - 1]) {
                std::string mutated = sequence;
                mutated[k - 1] = bp;
                MutatedKmerInterval snv(sequence, mutated);
                mutseqs.push_back(snv);
            }
        }
        assert(mutseqs.size() == 3);
    }

    void check_mutated_kmers(CountingHash& countgraph, bool debug)
    {
        for (auto& snv : mutseqs) {
            snv.check_mutated_kmers(countgraph, debug);
            assert(snv.kmer_abund.size() == k);
        }
    }
};

std::ostream& operator<<(std::ostream& os, const KmerInterval& ki)
{
    os << '>' << ki.sequence << '\n';
    for (int i = 0; i < ki.k; i++) {
        os << ' ';
    }
    os << '|' << '\n';

    for (auto& snv : ki.mutseqs) {
        bool first = true;
        for (auto kmerfreq : snv.kmer_abund) {
            if (first) {
                first = false;
            }
            else {
                os << ' ';
            }
            os << kmerfreq;
        }
        os << '\n';
    }

    return os;
}

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

std::vector<HashIntoType> get_n_primes_near_x(int x, int n)
{
    std::vector<HashIntoType> primes;
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

void count_mutation_collisions(CountingHash& countgraph, IParser *parser, int k,
                               unsigned long limit, bool debug)
{
    Sequence seq;
    unsigned long nucl_count = 0;

    while (!parser->is_complete()) {
        try {
            seq = parser->get_next_read();
        } catch (NoMoreReadsAvailable &e) {
            break;
        }

        for (unsigned long i = k - 1; i + k <= seq.sequence.length(); i++) {
            unsigned long minpos = i - k + 1;
            int subseqlen = (2*k) - 1;
            std::string interval = seq.sequence.substr(minpos, subseqlen);
            KmerInterval ki(interval, k);
            ki.check_mutated_kmers(countgraph, debug);
            std::cout << ki;

            nucl_count++;
            if (nucl_count % 100000 == 0) {
                std::cerr << "  processed " << (float)nucl_count / (float)1000
                          << " Kb of sequence" << std::endl;
            }
            if (limit > 0 && nucl_count > limit) {
                return;
            }
        }
    }
}

void print_usage(std::ostream& stream = std::cerr)
{
    stream << "Usage: snp-power [options] seqs.fa [genome.fa]\n";
    stream << "  options:\n";
    stream << "    -h    print this help message and exit\n";
    stream << "    -d    print debugging output\n";
    stream << "    -k    k-mer length (default: 31)\n";
    stream << "    -x    approx table size (default: 5e8)\n";
    stream << "    -N    num tables (default: 4)\n";
    stream << "    -l    limit number of positions to process\n";
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

    bool debug = false;
    int k = 31;
    int targetsize = 5e8;
    int numtables = 4;
    unsigned long limit = 0;

    char c;
    while ((c = getopt (argc, (char *const *)argv, "hdk:x:N:l:")) != -1) {
        if (c == 'h') {
            print_usage(std::cout);
            return 0;
        }
        else if (c == 'd') {
            debug = true;
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
            std::cerr << "Unknown option '" << c << "'\n";
            print_usage();
            return 1;
        }
    }

    std::string infile(argv[optind]);
    std::string refrfile(argv[optind]);
    if (argc > optind + 1) {
        refrfile = argv[optind + 1];
    }


    std::cerr << "allocating countgraph" << std::endl;
    std::vector<HashIntoType> tablesizes = get_n_primes_near_x(targetsize, numtables);
    CountingHash countgraph(k, tablesizes);

    std::cerr << "consuming input" << std::endl;
    unsigned int seqs_consumed = 0;
    unsigned long long kmers_consumed = 0;
    countgraph.consume_fasta(refrfile, seqs_consumed, kmers_consumed);
    std::cerr << "consumed " << seqs_consumed << " sequence(s) and "
              << kmers_consumed << " " << k << "-mers" << std::endl;

    std::cerr << "generating mutations" << std::endl;
    IParser *parser = IParser::get_parser(infile);
    count_mutation_collisions(countgraph, parser, k, limit, debug);
    delete parser;

    return 0;
}
