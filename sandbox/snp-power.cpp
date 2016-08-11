#include <assert.h>
#include <cmath>
#include <iostream>
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

void get_windows(std::string& sequence, int k, vstring& windows) {
    for (unsigned long i = 0; i < sequence.length(); i++) {
        if (i < k - 1 || i + k > sequence.length()) {
            continue;
        }
        int minpos = i - k + 1;
        int len = (2*k) - 1;
        std::string window = sequence.substr(minpos, len);
        windows.push_back(window);
    }
}

void get_mutations(Read& seq, int k, unsigned long limit, vstring& mutseqs)
{
    int i = k - 1;
    std::string nucleotides("ACGT");
    vstring windows;
    get_windows(seq.sequence, k, windows);

    int counter = 0;
    vstring::iterator it;
    for (it = windows.begin(); it < windows.end(); it++) {
        counter++;
        std::string window = *it;
        assert(window.length() == (2*k) - 1);

        if (counter > 0 && counter % 10000 == 0) {
            std::cerr << "Processed " << counter / 1000 << " Kb of sequence '"
                      << seq.name << "'" << std::endl;
        }
        if (limit > 0 && counter >= limit) {
            break;
        }

        for (auto nucl : nucleotides) {
            if (window[i] != nucl) {
                std::string mutseq = window;
                mutseq[i] = nucl;
                mutseqs.push_back(mutseq);
            }
        }
    }
}

void count_mutation_collisions(Hashbits& nodegraph, IParser *parser, int k, unsigned long limit)
{
    unsigned long total = 0;
    unsigned long hits = 0;
    Read seq;

    while (!parser->is_complete()) {
        try {
            seq = parser->get_next_read();
        } catch (NoMoreReadsAvailable &e) {
            break;
        }

        vstring mutseqs;
        get_mutations(seq, k, limit, mutseqs);
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
    std::cout << total << " " << hits << " " << (float)hits / (float)total << std::endl;
}

int main(int arg, const char **argv)
{
    assert(isprime(37));
    assert(!isprime(51));
    assert(isprime(15485863));
    assert(!isprime(15485865));

    int k = 31;
    int targetsize = 5e8;
    int numtables = 4;
    unsigned long limit = 10000;
    std::string infile(argv[1]);

    std::cerr << "allocating nodegraph" << std::endl;
    vint tablesizes = get_n_primes_near_x(targetsize, numtables);
    Hashbits nodegraph(k, tablesizes);

    std::cerr << "consuming input" << std::endl;
    unsigned int reads_consumed = 0;
    unsigned long long kmers_consumed = 0;
    nodegraph.consume_fasta(infile, reads_consumed, kmers_consumed);
    std::cerr << "consumed " << reads_consumed << " reads and "
              << kmers_consumed << " " << k << "-mers" << std::endl;

    std::cerr << "generating mutations" << std::endl;
    IParser *parser = IParser::get_parser(infile);
    count_mutation_collisions(nodegraph, parser, k, limit);
    delete parser;

    return 0;
}
