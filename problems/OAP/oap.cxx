#include <algorithm>
#include <iostream>
#include <fstream>
#include <utility>
#include <string>
#include <vector>
#include <tuple>
#include <cstdlib>

auto fasta(std::istream& in)
{
    std::vector< std::pair <std::string, std::string> > data;
    
    std::string line, label, sequence;
    bool started = false;
    while( std::getline(in, line).good() ) {
        if (line[0] == '>') {
            if (started)
                data.push_back( std::make_pair(label, sequence) );
            started = true;
            label = line.substr(1, line.length());
            sequence = "";
        } else {
            sequence += line;
        }
    }
    if (started)
        data.push_back( std::make_pair(label, sequence) );
    
    return data;
}

auto overlap_alignment(std::string seq1, std::string seq2, int match, int mismatch, int gap)
{
    enum move { start, ins, del, diag };
    const unsigned int n = seq1.length();
    const unsigned int m = seq2.length();
    
    std::vector< std::vector<move> > moves(n+1, std::vector<move>(m+1));
    std::vector< std::vector<int> > score(n+1, std::vector<int>(m+1));
    
    for (int i=0; i<=n; i++)
        for (int j=0; j<=m; j++) {
            moves[i][j] = start;
            score[i][j] = 0;
        }
    
    for (int i=1; i<=n; i++) {
        score[i][0] = i*gap;
        moves[i][0] = ins;
    }
    
    for (int i=1; i<=n; i++) {
        const char basei = seq1[i-1];
        for (int j=1; j<=m; j++) {
            int match_score = score[i-1][j-1];
            if (basei == seq2[j-1])
                match_score += match;
            else
                match_score += mismatch;
            
            const int ins_score = score[i-1][j] + gap;
            const int del_score = score[i][j-1] + gap;
            
            auto score_and_match = std::max(std::make_pair(ins_score, ins), std::make_pair(del_score, del));
            score_and_match = std::max(score_and_match, std::make_pair(match_score, diag));
            
            score[i][j] = score_and_match.first;
            moves[i][j] = score_and_match.second;
        }
    }

    // backtrack
    int maxscore = std::numeric_limits<int>::min();

    // highest i that has maximum score
    int maxloc = -1;
    for (int i=0; i<=n; i++) {
        if (score[i][m] >= maxscore) {
            maxloc = i;
            maxscore = score[i][m];
        }
    }
    
    int i = maxloc;
    int j = m;
    std::string aug1 = "";
    std::string aug2 = "";
    
    while (true) {
        if ((i<=0) || (j<=0))
            break;
        
        if (moves[i][j] == start)
            break;
        
        if (moves[i][j] == diag) {
            aug1 = seq1[i-1] + aug1;
            aug2 = seq2[j-1] + aug2;
            i--;
            j--;
        } else if (moves[i][j] == del) {
            aug1 = "-" + aug1;
            aug2 = seq2[j-1] + aug2;
            j--;
        } else {
            aug1 = seq1[i-1] + aug1;
            aug2 = "-" + aug2;
            i--;
        }
    }
    
    return std::tuple<int, std::string, std::string>(maxscore, aug1, aug2);
}

int main(int argc, char** argv)
{
    std::vector< std::pair <std::string, std::string> > fasta_data;
    if ( argc > 1 ) {
        std::ifstream fasta_in(argv[1]);
        fasta_data = fasta(fasta_in);
    } else {
        fasta_data = fasta(std::cin);
    }
    
    int gap_score = -2;
    int match_score = 1;
    int mismatch_score = -2;
    
    if (argc > 2)
        gap_score = atoi(argv[2]);
    if (argc > 3)
        match_score = atoi(argv[3]);
    if (argc > 4)
        mismatch_score = atoi(argv[4]);
    
    auto result = overlap_alignment(fasta_data[1].second, fasta_data[0].second,
                                    match_score, mismatch_score, gap_score);
    
    std::cout << std::get<0>(result) << std::endl;
    std::cout << std::get<2>(result) << std::endl;
    std::cout << std::get<1>(result) << std::endl;
}

