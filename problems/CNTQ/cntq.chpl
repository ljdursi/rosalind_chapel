module CNTQ {
  use BigInteger;
  use List;
  use Sort;
  use Search;
  const defaultfilename="none";
  config const infile=defaultfilename;
  config const modulus=1000000;

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
        tokens.push_back(new Token(TokenType.openp, ""));
        continue;
      }
     
      if item(1) == ")" {
        var name = item(2..);
        tokens.push_back(new Token(TokenType.closep, name));
        continue;
      }

      var name = item;
      tokens.push_back(new Token(TokenType.node, name));
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
    var namedprogeny: domain(string);

    proc init(tokens: [] Token) {
      const low = tokens.domain.dim(1).low;
      const high = tokens.domain.dim(1).high;
      const n = tokens.size;

      this.complete();
      if (n > 0) {
        this.name = tokens[high].name;
      }

      if (n > 1) {
        const sublists = partition_list(tokens[low+1..high-1]);
        for sublist in sublists {
          const subtree = new unmanaged Tree(sublist.arr);
          this.children.push_back(subtree);

          for name in subtree.namedprogeny do
            this.namedprogeny.add(name);
        }
      }
      if this.name != "" then
        this.namedprogeny.add(name);
    }

    proc nchildren(): int {
      return this.children.size;
    }
  
    iter edgewalks() : 2*Tree {
      for child in children {
        yield (this, child);
        for edge in child.edgewalks() do
          yield edge;
      }
    }

    proc edgesets(child, allleafs) {
      var childsets : [1..2] domain(string);
      var othersets : [1..2] domain(string);
      var otherchildren : [1..0] Tree;

      for c in this.children do
        if c != child then
          otherchildren.push_back(c);

      for i in 1..child.nchildren() do
        childsets[i] = child.children[i].namedprogeny;

      othersets[1] = otherchildren[1].namedprogeny;
      if otherchildren.size > 1 then
        othersets[2] = otherchildren[2].namedprogeny;
      else
        othersets[2] = allleafs - othersets[1] - childsets[1] - childsets[2];

      return (childsets[1], childsets[2], othersets[1], othersets[2]);
    }

    proc nquartets() {
      const allleafs = this.namedprogeny;
      var nq = new bigint(0);

      for (edgep, edgec) in this.edgewalks() {
        const (c1, c2, o1, o2) = edgep.edgesets(edgec, allleafs);
        const csize = c1.size + c2.size;
        const osize = o1.size + o2.size;

        var nlocal = c1.size*c2.size*(osize*(osize-1))/2;
        nlocal += o1.size*o2.size*(csize*(csize-1))/2;
        nq += nlocal;
      }

      nq = nq / 2;
      return (nq % modulus) : int;
    }
  }
  
  proc main() {
    try! {
      var lines: [1..0] string;
      for line in lines_from_file(infile, defaultfilename) do
        lines.push_back(line);

      const tokens = tokenize(lines[2]);
      const phylogeny = new owned Tree(tokens);

      const nq = phylogeny.nquartets();
      writeln(nq);
    }
  }
}
