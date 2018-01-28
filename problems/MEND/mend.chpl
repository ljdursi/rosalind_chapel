module MEND {
  use List;
  use Sort;
  use Search;
  use LinearAlgebra;
  const defaultfilename="none";
  config const infile=defaultfilename;

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

  proc allele_count_probabilities(parent1, parent2) {
    var T0 = Matrix(3, 3);  
    var T2 = Matrix(3, 3);  
    var allele_probs = Vector(3), curprobs = Vector(3);

    for i in 0..2 do
      for j in 0..2 {
        T0[i, j] = (2-i)*(2-j)/4.0;
        T2[i, j] = i*j/4.0;
      }

    allele_probs[0] = dot( parent2, dot( T0, parent1 ));
    allele_probs[2] = dot( parent2, dot( T2, parent1 ));
    allele_probs[1] = 1.0 - allele_probs[0] - allele_probs[2];
    return allele_probs;
  }

  proc genotype_to_probabilities(gt) {
    const prob_dict = ['AA' => 0, 'Aa' => 1, 'aA' => 1, 'aa' => 2];

    var p = Vector(3);
    p[prob_dict[gt]] = 1.0;

    return p;
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

  class Tree {
    var name: string;
    var children: [1..0] Tree;
    var gt_probs = Vector(3);

    proc init(tokens: [] Token) {
      const low = tokens.domain.dim(1).low;
      const high = tokens.domain.dim(1).high;
      const n = tokens.size;

      super.init();
      if (n > 0) {
        this.name = tokens[high].name;
        if this.name != "" then
          this.gt_probs = genotype_to_probabilities(this.name);
      }

      if (n > 1) {
        const sublists = partition_list(tokens[low+1..high-1]);
        for sublist in sublists {
          const subtree = new Tree(sublist.arr);
          this.children.push_back(subtree);
        }
      }

      if this.name == "" then
        if this.children.size == 2 then
          this.gt_probs = allele_count_probabilities(this.children[1].gt_probs,
                                                     this.children[2].gt_probs);
    }
  
    iter these() : Tree {
      for child in children do
        for e in child.these() do
          yield e;
      yield this;
    }
  }

  proc main() {
    try! {
      for line in lines_from_file(infile, defaultfilename) {
        const tokens = tokenize(line);
        const phylogeny = new Tree(tokens);
        
        writeln(phylogeny.gt_probs);

        delete phylogeny;
      }
    }
  }
}
