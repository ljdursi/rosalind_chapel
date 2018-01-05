module NWCK {
  use Regexp;
  const defaultfilename="none";
  config const infile=defaultfilename;

  enum Token { openp, closep, close_w_parent, done, comma, node };

  proc strings_from_file(filename: string, defaultfilename: string) throws {
    var line: string;
    var lines: [1..0] string;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    while (channel.readline(line)) {
      lines.push_back(line.strip());
    }

    return lines;
  }

  proc tokenize(newick) throws {
    var parser = compile("(;|,|\\(|\\)[^\\(\\),;]+|\\)|[^\\(\\),;]+)");
    var tokens : [1..0] (Token, string);

    const fixed_tokens = {',', ';', '(', ')'};
    const token_fm_string = [',' => Token.comma, ';' => Token.done, 
                             '(' => Token.openp, ')' => Token.closep];

    for item in parser.split(newick) {
      if item == "" then
        continue;
      if fixed_tokens.member(item) {
        tokens.push_back((token_fm_string[item], item));
        continue;
      }
      if item(1) == ")" then
        tokens.push_back((Token.close_w_parent, item(2..)));
      else
        tokens.push_back((Token.node, item));
    }
    return tokens;
  }

  proc find_dist(newick, node1, node2) throws {
    const tokens = tokenize(newick);
    var start = node1;
    var finish = node2;
    var loc1, loc2 : int;

    for i in 1..tokens.size {
      if tokens[i][2] == node1 then
        loc1 = i;
      if tokens[i][2] == node2 then
        loc2 = i;
    }

    if loc2 < loc1 {
        (loc2, loc1) = (loc1, loc2);
        (start, finish) = (node2, node1);
    }

    var state : [1..0] Token;
    state.push_back(Token.node);
    
    for (token_type, token_val) in tokens[loc1+1..loc2] {
      if token_type == Token.closep || token_type == Token.close_w_parent {
        if state[state.size] == Token.openp then
          state.pop_back();
        else 
          state.push_back(token_type);
      }
      if token_type == Token.openp then
        state.push_back(token_type);
    }

    if tokens[loc2][1] == Token.close_w_parent then
      return state.size - 1;
    return state.size + 1;
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
