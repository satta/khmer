/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2016, The Regents of the University of California.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Michigan State University nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
LICENSE (END)

Contact: khmer-project@idyll.org
*/
#include <seqan/seq_io.h> // IWYU pragma: keep
#include <seqan/sequence.h> // IWYU pragma: keep
#include <seqan/stream.h> // IWYU pragma: keep
#include <fstream>

#include "khmer_exception.hh"
#include "seqio.hh"

namespace khmer
{

namespace sequence_io
{

// -----------------------------------------------------------------------------
// Parser base class
// -----------------------------------------------------------------------------

Parser::Parser() : _num_reads(0)
{
    int haderror;

    haderror = regcomp(
        &_read_1_pattern,
        "^.+(/1| 1:[YN]:[[:digit:]]+:[[:alpha:]]+).{0}",
        REG_EXTENDED
    );
    if (haderror) {
        throw khmer_exception("Could not compile R1 regex");
    }

    haderror = regcomp(
        &_read_2_pattern,
        "^.+(/2| 2:[YN]:[[:digit:]]+:[[:alpha:]]+).{0}",
        REG_EXTENDED
    );
    if (haderror) {
        throw khmer_exception("Could not compile R2 regex");
    }
}

Parser::~Parser()
{
    regfree(&_read_1_pattern);
    regfree(&_read_2_pattern);
}

void Parser::imprint_next_read_pair(ReadPair& pair, PairMode mode)
{
    if (mode == PAIR_MODE_IGNORE_UNPAIRED) {
        _imprint_next_read_pair_in_ignore_mode(pair);
    } else if (mode == PAIR_MODE_ERROR_ON_UNPAIRED) {
        _imprint_next_read_pair_in_error_mode(pair);
    }
#if (0)
    else if (mode == PAIR_MODE_ERROR_ON_UNPAIRED) {
        _imprint_next_read_pair_in_allow_mode(pair);
    }
#endif
    else {
        std::string msg = "Unknown pair reading mode: " + std::to_string(mode);
        throw UnknownPairReadingMode(msg);
    }
}

#if (0)
void Parser::_imprint_next_read_pair_in_allow_mode(ReadPair& pair)
{
    // TODO: Implement.
    //	     Probably need caching of reads between invocations
    //	     and the ability to return pairs which are half empty.
}
#endif

void Parser::_imprint_next_read_pair_in_ignore_mode(ReadPair& pair)
{
    Read& read_1 = pair.first;
    Read& read_2 = pair.second;
    regmatch_t match_1, match_2;

    // Hunt for a read pair until one is found or end of reads is reached.
    while (true) {

        // Toss out all reads which are not marked as first of a pair.
        // Note: We let any exception, which flies out of the following,
        //	 pass through unhandled.
        while (true) {
            imprint_next_read(read_1);
            int result = regexec(&_read_1_pattern, read_1.name.c_str(), 1,
                &match_1, 0
            );
            if (result == REG_NOMATCH) {
                break;
            }
        }

        // If first read of a pair was found, then insist upon second read.
        // If not found, then restart search for pair.
        // If found, then validate match.
        // If invalid pair, then restart search for pair.
        imprint_next_read(read_2);
        int result = regexec(&_read_2_pattern, read_2.name.c_str(), 1,
            &match_2, 0
        );
        if (result == REG_NOMATCH) {
            if (_is_valid_read_pair(pair, match_1, match_2)) {
                break;
            }
        }

    }
}

void Parser::_imprint_next_read_pair_in_error_mode(ReadPair& pair)
{
    Read& read_1 = pair.first;
    Read& read_2 = pair.second;
    regmatch_t match_1, match_2;

    // Note: We let any exception, which flies out of the following,
    //	     pass through unhandled.
    imprint_next_read(read_1);
    imprint_next_read(read_2);

    int r1 = regexec(&_read_1_pattern, read_1.name.c_str(), 1, &match_1, 0);
    int r2 = regexec(&_read_2_pattern, read_2.name.c_str(), 1, &match_2, 0);
    if (r1 == REG_NOMATCH || r2 == REG_NOMATCH) {
        throw InvalidReadPair();
    }

    if (!_is_valid_read_pair(pair, match_1, match_2)) {
        throw InvalidReadPair();
    }
}

bool Parser::_is_valid_read_pair(ReadPair& pair, regmatch_t& m1, regmatch_t& m2)
{
    bool smatch = m1.rm_so == m2.rm_so;
    bool ematch = m1.rm_eo == m2.rm_eo;
    bool submatch = pair.first.name.substr(0, m1.rm_so) ==
                    pair.second.name.substr(0, m2.rm_so);
    return smatch && ematch && submatch;
}

// -----------------------------------------------------------------------------
// Fasta/Fastq parser
// -----------------------------------------------------------------------------

struct FastxParser::Handle
{
    seqan::SequenceStream stream;
    uint32_t seqan_spin_lock;
};

FastxParser::FastxParser(char const * filename) : Parser()
{
    _handle = new FastxParser::Handle();
    seqan::open(_handle->stream, filename);
    if (!seqan::isGood(_handle->stream)) {
        std::string message = "Could not open ";
        message = message + filename + " for reading.";
        throw InvalidStream(message);
    } else if (seqan::atEnd(_handle->stream)) {
        std::string message = "File ";
        message = message + filename + " does not contain any sequences!";
        throw InvalidStream(message);
    }
    __asm__ __volatile__ ("" ::: "memory");
    _handle->seqan_spin_lock = 0;
}

FastxParser::~FastxParser()
{
    seqan::close(_handle->stream);
    delete _handle;
}

bool FastxParser::is_complete()
{
    return !seqan::isGood(_handle->stream) || seqan::atEnd(_handle->stream);
}

void FastxParser::imprint_next_read(Read& read)
{
    read.reset();
    int ret = -1;
    const char *invalid_read_exc = NULL;
    while (!__sync_bool_compare_and_swap(& _handle->seqan_spin_lock, 0, 1));
    bool atEnd = seqan::atEnd(_handle->stream);
    if (!atEnd) {
        ret = seqan::readRecord(read.name, read.sequence, read.quality,
                                _handle->stream);
        if (ret == 0) {
            // Detect if we're parsing something w/ qualities on the first read
            // only
            if (_num_reads == 0 && read.quality.length() != 0) {
                _have_qualities = true;
            }

            // Handle error cases, or increment number of reads on success
            if (read.sequence.length() == 0) {
                invalid_read_exc = "Sequence is empty";
            } else if (_have_qualities && (read.sequence.length() != \
                                           read.quality.length())) {
                invalid_read_exc = "Sequence and quality lengths differ";
            } else {
                _num_reads++;
            }
        }
    }
    __asm__ __volatile__ ("" ::: "memory");
    _handle->seqan_spin_lock = 0;
    // Throw any error in the read, even if we're at the end
    if (invalid_read_exc != NULL) {
        throw InvalidRead(invalid_read_exc);
    }
    // Throw NoMoreReadsAvailable if none of the above errors were raised, even
    // if ret == 0
    if (atEnd) {
        throw NoMoreReadsAvailable();
    }
    // Catch-all error in readRecord that isn't one of the above
    if (ret != 0) {
        throw StreamReadError();
    }
}

} // namespace read_parsers

} // namespace khmer
