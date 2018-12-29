module LREP {
  use List;
  const defaultfilename="none";
  config const infile=defaultfilename;

  // use stdin if filename == defaultfilename
  proc params_from_file(filename: string, defaultfilename: string) throws {
    var channel = stdin;

    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var fullstring: string;
    var minrepeat: uint;
    var success = channel.readline(fullstring);
    fullstring.strip();

    success = channel.readln(minrepeat);

    var lines: list(string);
    var line: string;
    var firstline: string;

    while (channel.readline(line)) {
      if firstline == "" then
        firstline = line;
      lines.push_back(line.strip());
    }

    const items = firstline.split();
    return (fullstring, minrepeat, items[1], lines);
  }

  class Tree {
    var start: uint;
    var length: uint;
    var totallength: uint;
    var nendings: uint;
    var name: string;
    var children: [1..0] unmanaged Tree;
  
    iter these() : Tree {
      for child in children do
        for e in child.these() do
          yield e;
      yield this;
    }

    proc deinit() {
      for child in children do
        delete child;
    }

  }
  
  proc partition_list(all_lines: list(string), children: [?D] string) {
    const nchildren = children.size : uint;
    var sublists: [1..#nchildren] list(string);
  
    var nodes: domain(string);
    var subtrees: [nodes] uint;
  
    for i in 1..#nchildren {
      nodes += children[i:int];
      subtrees[children[i:int]] = i:uint;
    }
  
    for line in all_lines {
      const items = line.split();
      const parent_name = items[1];
      const child_name = items[2];
      if nodes.member(parent_name) {
        const subtree_num = subtrees[parent_name];
        sublists[subtree_num:int].append(line);
        subtrees[child_name] = subtree_num;
      }
    }
  
    return sublists;
  }
  
  proc buildTree(ref lines: list(string), fullstring: string,
                 name: string, start: uint = 0,
                 length: uint = 0, totallength: uint = 0) : unmanaged Tree throws {
    const l = fullstring.length;
    var children_names: [1..0] string;
    var children_starts: [1..0] uint;
    var children_lens: [1..0] uint;
  
    var tree = new unmanaged Tree(start, length, totallength, 0: uint, name);
    tree.nendings = 0;
    if start + length >= l then
      tree.nendings = 1;

    if lines.size == 0 then
      return tree;
  
    var line = lines.pop_front();
    var items = line.split(); 
    while (items[1] == name) {
      try! {
        const parent_name = items[1];
        const child_name = items[2];
        const start = items[3] : uint;
        const len = items[4] : uint;
    
        children_names.push_back(child_name);
        children_starts.push_back(start);
        children_lens.push_back(len);
    
        if lines.size > 0 {
          line = lines.pop_front();
          items = line.split(); 
        } else 
          break;
      }
    } 
    lines.push_front(line);
  
    var sublists = partition_list(lines, children_names);
    for (sublist, subname, substart, sublen) in zip(sublists, children_names, children_starts, children_lens) {
        var child = buildTree(sublist, fullstring, subname, substart, sublen, sublen + totallength);
        tree.nendings += child.nendings;
        tree.children.push_back(child);
    }
    return tree;
  }
  


  proc main() {
    try! {
      var (totstring, minrepeats, rootname, lines) = params_from_file(infile, defaultfilename);
      const suffixtree = buildTree(lines, totstring, rootname);

      var maxlen = 0 :uint;
      var maxstart = 0 :uint;
      var maxname = "";

      for subtree in suffixtree {
        if subtree.nendings >= minrepeats {
          if subtree.totallength > maxlen {
            maxlen = subtree.totallength;
            maxstart = subtree.start - subtree.totallength + subtree.length;
            maxname = subtree.name;
          }
        }
      }
      writeln(totstring(maxstart..#maxlen));
      delete suffixtree;
    }
  }
}
