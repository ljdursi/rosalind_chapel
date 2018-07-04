module NKEW {
  use Regexp;
  const defaultfilename="none";
  config const infile=defaultfilename;

  enum TokenType { openp, closep, done, comma, node };

  record Token {
    var tt: TokenType;
    var name: string;
    var weight: int;
  }
  
  iter strings_from_file(filename: string, defaultfilename: string) throws {
    var line: string;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    while (channel.readline(line)) do
      yield line.strip();
  }

  proc parseweight(name: string) : (string, int) throws {
    if name.length == 0 then
      return ("", 0);

    var items = name.split(':');
    if items.size == 1 then
      return (items[1], 1);
    else {
      try! {
        const wt: int = items[2]:int;
        return (items[1], wt);
      }
    }
  }

  proc tokenize(newick) throws {
    var parser = compile('(;|,|\\(|\\)[^\\(\\),;]*|[^\\(\\),;]+)');
    var tokens : [1..0] Token;
    var last = -1;
    var i = 0;

    const fixed_tokens = {',', ';', '('};
    const token_fm_string = [',' => TokenType.comma, ';' => TokenType.done, 
                             '(' => TokenType.openp];

    for item in parser.split(newick) {
      if item == "" then
        continue;

      i += 1;
      if fixed_tokens.member(item) {
        tokens.push_back(new Token(token_fm_string[item], "", 0));
        if item == ";" then
          last = i;
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

    return tokens[1..last-1];
  }

  proc find_dist(newick, node1, node2) throws {
    const tokens = tokenize(newick);
    var start = node1;
    var finish = node2;
    var loc1, loc2 : int;

    for i in 1..tokens.size {
      if tokens[i].name == node1 then
        loc1 = i;
      if tokens[i].name == node2 then
        loc2 = i;
    }

    if loc2 < loc1 {
        (loc2, loc1) = (loc1, loc2);
        (start, finish) = (node2, node1);
    }

    var path : [1..0] Token;
    path.push_back(tokens(loc1));

    var level = 0;
    var minlevel = 0;
    
    for (idx, token) in zip((loc1+1).., tokens[loc1+1..]) {
      const tt = token.tt;
      const name = token.name;
      const weight = token.weight;

      if tt == TokenType.closep {
        level -= 1;
        if path[path.size].tt == TokenType.openp then
          path.pop_back();
        else 
          path.push_back(token);
      } else if tt == TokenType.openp {
        level += 1;
        path.push_back(token);
      } else if name == start || name == finish {
        path.push_back(token);
      }

      // we may have to keep seeing )s to find remaining edge weights
      if level < minlevel && idx < loc2 then
        minlevel = level;
      if idx >= loc2 && level <= minlevel then
        break;
    }

    var pathdist = 0;
    for token in path do
      pathdist += token.weight;

    if path[path.size].name == node2 && level < minlevel then
        pathdist -= path[path.size].weight;

    return pathdist;
  }


  proc main() {
    try! {
      var lines = strings_from_file(infile, defaultfilename);
      var distances : [1..0] int;
      for (newick, line2) in zip(lines[1.. by 3], lines[2.. by 3]) {
        var items = line2.split();

        distances.push_back(find_dist(newick, items[1], items[2]));
      }
      writeln(distances);
    }
  }
}
