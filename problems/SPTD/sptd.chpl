module SPTD {
  use List;
  use Sort;
  use Search;
  const defaultfilename="none";
  config const infile=defaultfilename;

  // use stdin if filename == defaultfilename
  proc lines_from_file(filename: string, defaultfilename: string) throws {
    var channel = stdin;

    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var line: string;
    var lines: [1..0] string;
    while (channel.readline(line)) do
      lines.push_back(line.strip());
  
    return lines;
  }

  enum TokenType { openp, closep, node };

  record Token {
    var tt: TokenType;
    var name: string;
  }

  proc tokenize(newick) throws {
    var parser = compile('(;|,|\\(|\\)[^\\(\\),;]*|[^\\(\\),;]+)');
    var tokens : [1..0] Token;

    for item in parser.split(newick) {
      if item == "" then
        continue;

      if item == "," then
        continue;

      if item == ";" then
        break;

      if item == "(" {
        tokens.push_back(new Token(TokenType.openp, ""));
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
    var namedprogeny: [1..0] string;

    proc init(tokens: [] Token) {
      const low = tokens.domain.dim(1).low;
      const high = tokens.domain.dim(1).high;
      const n = tokens.size;

      super.init();
      if (n > 0) {
        this.name = tokens[high].name;
      }

      if (n > 1) {
        const sublists = partition_list(tokens[low+1..high-1]);
        for sublist in sublists {
          const subtree = new Tree(sublist.arr);
          this.children.push_back(subtree);
          if subtree.name != "" then
            this.namedprogeny.push_back(subtree.name);

          for name in subtree.namedprogeny do
            this.namedprogeny.push_back(name);
        }
      }
    }
  
    iter these() : Tree {
      for child in children do
        for e in child.these() do
          yield e;
      yield this;
    }

    proc splits(allnames) {
      var splits_set : domain(string);

      proc as_string(d: domain(string)) {
        var result = "";
        for item in d.sorted() do
          result += item + " ";
        return result;
      }

      const n = allnames.size;
      const allnames_set : domain(string) = allnames;

      for subtree in this {
        var subnames_set : domain(string) = subtree.namedprogeny;
        const subn = subnames_set.size;
        if subn <= 1 || subn >= n-1 then
          continue;
        if subnames_set.size > n - subnames_set.size then
          subnames_set = allnames_set - subnames_set;

        splits_set += as_string(subnames_set);
      }
      return splits_set;
    }
  }

  
  proc main() {
    try! {
      const lines = lines_from_file(infile, defaultfilename);
      const allnames = lines[1].split();

      const phylogeny1 = new Tree(tokenize(lines[2]));
      const phylogeny2 = new Tree(tokenize(lines[3]));

      const diff = phylogeny1.splits(allnames) ^ phylogeny2.splits(allnames);
      writeln(diff.size);
    }
  }
}
