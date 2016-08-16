// g++ --std=c++11 -O3 -I ../lib -I ../third-party/seqan/core/include/ -o snp-power snp-power.cpp ../lib/liboxli.a

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
    std::string& mutseq;
    unsigned long unique_kmer_count;
    vint kmer_abund;

    MutatedKmerInterval(std::string& sequence, std::string& mutatedseq)
        : original(sequence), mutseq(mutatedseq) {}

    unsigned long total_kmer_count() {
        return kmer_abund.size();
    }

    unsigned long non_unique_kmer_count() {
        return total_kmer_count() - unique_kmer_count;
    }

    void check_mutated_kmers(CountingHash& countgraph, bool debug)
    {
        vstring kmers;
        countgraph.get_kmers(mutseq, kmers);
        int kmer_offset = 0;
        for (auto kmer : kmers) {
            int kmerfreq = countgraph.get_count(kmer.c_str());
            kmer_abund.push_back(kmerfreq);
            if (kmerfreq == 0) {
                unique_kmer_count++;
            }
            else {
                if (debug) {
                    std::cerr << "DEBUG orig seq: " << original << std::endl;
                    std::cerr << "DEBUG mutd seq: " << mutseq << std::endl;
                    std::cerr << "DEBUG         : ";
                    for (int i = 0; i < kmer_offset; i++) {
                        std::cerr << ' ';
                    }
                    std::cerr << kmer << std::endl << std::endl;
                }
            }
            kmer_offset++;
        }
    }
};
typedef std::vector<MutatedKmerInterval> vmut;

/// A nucleotide of interest, along with k - 1 nucleotides on either side; in
/// other words, the interval containing all k-mers overlapping with the
/// nucleotide of interest.
class KmerInterval
{
public:
    Sequence& seq;
    int k;
    int i;
    std::string nucl;

    KmerInterval(Sequence& sequence, int ksize)
        : seq(sequence), k(ksize), i(k - 1), nucl("ACGT") {}

    bool next_mutated_seqs(vmut& mutseqs)
    {
        mutseqs.clear();
        if (i + k > seq.sequence.length()) {
            return false;
        }

        int minpos = i - k + 1;
        int len = (2*k) - 1;
        std::string subseq = seq.sequence.substr(minpos, len);
        for (auto bp : nucl) {
            std::string mutated = subseq;
            mutated[k - 1] = bp;
            MutatedKmerInterval snv(subseq, mutated);
            mutseqs.push_back(snv);
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

    int kmers_hit = 0;
    int kmers_unique = 0;
    std::vector<int> uniq_k_per_snv;
    std::vector<int> uniq_k_per_nucl;

    while (!parser->is_complete()) {
        try {
            seq = parser->get_next_read();
        } catch (NoMoreReadsAvailable &e) {
            break;
        }

        KmerInterval m(seq, k);
        vmut mutseqs;
        while (m.next_mutated_seqs(mutseqs) && (limit == 0 || nucl_count < limit)) {
            nucl_count++;
            if (nucl_count % 100000 == 0) {
                std::cerr << "  processed " << (float)nucl_count / (float)1000
                          << " Kb of sequence" << std::endl;
            }

            assert(mutseqs.size() == 3);
            int nucl_uniq_kmers = 0;
            for (auto mutseq : mutseqs) {
                mutseq.check_mutated_kmers(countgraph, debug);
                kmers_unique += mutseq.unique_kmer_count;
                nucl_uniq_kmers += mutseq.unique_kmer_count;
                kmers_hit += mutseq.non_unique_kmer_count();
                uniq_k_per_snv.push_back(mutseq.unique_kmer_count);
            }
            uniq_k_per_nucl.push_back(nucl_uniq_kmers);
        }
    }
    std::cout << (kmers_hit + kmers_unique) << " " << kmers_hit << " " << (float)kmers_hit / (float)(kmers_hit + kmers_unique)
              << std::endl;
}

void print_usage(std::ostream& stream = std::cerr)
{
    stream << "Usage: snp-power [options] seqs.fa [genome.fa]" << std::endl;
    stream << "  options:" << std::endl;
    stream << "    -h    print this help message and exit" << std::endl;
    stream << "    -d    print debugging output" << std::endl;
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
