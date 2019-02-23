module ALPH {
  // Alignment-Based Phylogeny
  // http://rosalind.info/problems/alph/
  use Regexp;
  use LinearAlgebra;
  const defaultfilename="none";
  config const infile=defaultfilename;

  const BASE_TO_IDX = ['A' => 1, 'C' => 2, 'G' => 3, 'T' => 4, '-' => 5];
  const IDX_TO_BASE = [1 => 'A', 2 => 'C', 3 => 'G', 4 => 'T', 5 => '-'];
  const NBASES = BASE_TO_IDX.domain.size;
  const PENALTY: [1..NBASES, 1..NBASES] int = 1 - eye(NBASES, int);

  // use stdin if filename == defaultfilename
  iter lines_from_file(filename: string, defaultfilename: string) throws {
    var channel = stdin;

    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var line: string;
    while (channel.readline(line)) do
      yield line.strip();
  }

  proc fasta_to_dict(lines) {
    // 
    // simple fasta parser
    //
    var started = false;
    var sequence = "", seqlabel = "";
    var labels: domain(string);
    var sequence_dict: [labels] string;

    for line in lines {
      if (line[1] == '>') {
        if (started) {
          labels += seqlabel;
          sequence_dict[seqlabel] = sequence;
        }
        started = true;
        seqlabel = line(2..).strip();
        sequence = "";
      } else {
        sequence += line.strip();
      }
    }
    if started {
      labels += seqlabel;
      sequence_dict[seqlabel] = sequence;
    }
    return sequence_dict;
  }

  enum TokenType { openp, closep, node };

  record Token {
    var tt: TokenType;
    var name: string;
    var weight: int;
  }

  proc tokenize(newick) throws {
    var parser = compile('(;|,|\\(|\\)[^\\(\\),;]*|[^\\(\\),;]+)');
    var tokens : [1..0] Token;
    var last = -1;
    var i = 0;

    for item in parser.split(newick) {
      if item == "" then
        continue;

      if item == "," then
        continue;

      if item == ";" then
        break;

      if item == "(" {
        tokens.push_back(new Token(TokenType.openp, "", 0));
        continue;
      }
     
      if item(1) == ")" {
        tokens.push_back(new Token(TokenType.closep, item(2..)));
        continue;
      }

      tokens.push_back(new Token(TokenType.node, item));
    }

    return tokens;
  }

  record tokenarray {
    var arr: [1..0] Token;
  }

  proc partition_list(tokens: [?] Token) {
    const low = tokens.domain.dim(1).low;
    const high = tokens.domain.dim(1).high;
    var dividers: [1..0] 2*int;
    var start = -1;
    var level = 0;

    for (i, token) in zip(low..high, tokens) {
      select (token.tt) {
        when (TokenType.openp) {
          level += 1;
          if level == 1 then
            start = i;
        }
        when (TokenType.closep) {
          level -= 1;
          if level == 0 {
            dividers.push_back((start, i));
            start = -1;
          }
        }
        when (TokenType.node) {
          if level == 0 then
            dividers.push_back((i, i));
        }
      }
    }

    var retval : [1..dividers.size] tokenarray;
    for (i, divider) in zip(1.., dividers) {
      const (a, b) = divider;
      for item in tokens[a..b] do
        retval[i].arr.push_back(item);
    }

    return retval;
  }

  proc hamming(seq1, seq2) {
    var dist = 0;
    for (base1, base2) in zip(seq1, seq2) do
      if base1 != base2 then
        dist += 1;

    return dist;
  }

  record assoc_array {
    var d: domain(string);
    var map: [d] string;
  }

  class Tree {
    const k: int;
    var name: string;
    var sequence: string;
    var children: [1..0] Tree;
    var seqmatrix_domain: domain(2);
    var seqmatrix: [seqmatrix_domain] int;

    proc init(tokens, sequences, k) {
      const low = tokens.domain.dim(1).low;
      const high = tokens.domain.dim(1).high;
      const n = tokens.size;
      super.init();

      this.k = k;
      this.seqmatrix_domain = {1..k, 1..NBASES};
      this.seqmatrix = 0;

      if (n > 0) then
        this.name = tokens[high].name;

      if sequences.domain.member(this.name) then
        this.sequence = sequences[this.name];

      if (n > 1) {
        const sublists = partition_list(tokens[low+1..high-1]);
        for sublist in sublists {
          var subtree = new unmanaged Tree(sublist.arr, sequences, k);
          this.children.push_back(subtree);
        }
      }
    }

    proc postinit() {
      if this.sequence != "" then
        this.seqmatrix[1..k, 1..NBASES] = this._sequence_to_matrix();
      else {
        for child in this.children do
          this.seqmatrix += this._project_matrix(child.seqmatrix);
      }
    }

    proc _sequence_to_matrix() {
      const l = sequence.length;
      var matrix : [1..l, 1..NBASES] int = 1;
      for (idx, base) in zip(1.., sequence) do
        matrix[idx, BASE_TO_IDX[base]] = 0;
      return matrix;
    }

    proc _project_matrix(mtx) {
      const ones : [1..NBASES] int = 1;
      var newmtx : [1..k, 1..NBASES] int = 0;
      for i in 1..k {
        const row = mtx[i, 1..NBASES];
        const choices : [PENALTY.domain] int = outer(ones, row) + PENALTY;
        for j in 1..NBASES do
          newmtx[i, j] = min reduce choices[j, 1..NBASES];
      }
      return newmtx;
    }

    proc _seqmatrix_to_sequence() {
      if sequence != "" then
        return;
      for i in 1..k {
        const (dmy, idx) = minloc reduce zip(seqmatrix[i, 1..NBASES], 1..NBASES);
        sequence += IDX_TO_BASE[idx];
      }
    }

    proc set_sequences() {
      _seqmatrix_to_sequence();
      const matrix = _sequence_to_matrix();
      for child in children {
        child.seqmatrix += matrix;
        child.set_sequences();
      }
    }

    proc calc_hamming() : int {
      var tree_hamming : int = 0;
      for child in children {
        tree_hamming += child.calc_hamming();
        tree_hamming += hamming(sequence, child.sequence);
      }

      return tree_hamming;
    }

    proc sequences() : assoc_array {
      var sequences_domain: domain(string);
      var sequences_dict: [sequences_domain] string;
      sequences_dict[name] = sequence;

      for child in children do {
        const child_aa = child.sequences();
        for lab in child_aa.d do
          sequences_dict[lab] = child_aa.map[lab];
      }

      return new assoc_array(sequences_domain, sequences_dict);
    }
  }

  proc main() {
    try! {
      const lines = lines_from_file(infile, defaultfilename);
      const tokens = tokenize(lines[1]);
      const given_sequences = fasta_to_dict(lines[2..]);

      var k = 0;
      for lab in given_sequences.domain do
        if k == 0 then
          k = given_sequences[lab].length;
        else
          assert(given_sequences[lab].length == k);

      var phylogeny = new owned Tree(tokens, given_sequences, k);
      phylogeny.set_sequences();

      const tree_sequences_aa = phylogeny.sequences();
      const score = phylogeny.calc_hamming();

      writeln(score);
      for taxon in tree_sequences_aa.d.sorted() do
        if !given_sequences.domain.member(taxon) {
          writeln(">", taxon);
          writeln(tree_sequences_aa.map[taxon]);
        }
    }
  }
}
