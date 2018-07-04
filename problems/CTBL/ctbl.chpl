module CTBL {
  use List;
  use Sort;
  use Search;
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

  enum TokenType { openp, closep, node };

  record Token {
    var tt: TokenType;
    var name: string;
    var weight: int;
  }

  proc parseweight(name: string) : (string, int) throws {
    if name.length == 0 then
      return ("", 0);

    var items = name.split(':');
    if items.size == 1 then
      return (items[1], 1);
    else
      try! {
        const wt = items[2]:int;
        return (items[1], wt);
      }
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
        var (name, weight) = parseweight(item(2..));
        tokens.push_back(new Token(TokenType.closep, name, weight));
        continue;
      }

      var (name, weight) = parseweight(item);
      tokens.push_back(new Token(TokenType.node, name, weight));
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
  }

  proc character_entry(allnames, subnames) {
    var result = "";
    var insym = "1", outsym = "0";
    if (allnames[1] == subnames[1]) {
      insym = "0";
      outsym = "1";
    }

    for name in allnames {
      const (found, dummy) = binarySearch(subnames, name);
      if found then
        result += insym;
      else
        result += outsym;
    }
    return result;
  }
  
  proc main() {
    try! {
      for line in lines_from_file(infile, defaultfilename) {
        const tokens = tokenize(line);
        const phylogeny = new Tree(tokens);
        
        const allnames = sorted(phylogeny.namedprogeny);
        const n = allnames.size;
        const trivial = {0, 1, n-1, n};
        for subtree in phylogeny {
          const subnames = sorted(subtree.namedprogeny);
          const nsub = subnames.size;

          if !trivial.member(nsub) then
            writeln(character_entry(allnames, subnames));
        }
      }
    }
  }
}
