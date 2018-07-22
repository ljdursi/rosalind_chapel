module MREP {
  const defaultfilename="none";
  config const infile=defaultfilename;
  config const threshold=20;

  record stringarray {
    var D: domain(1);
    var strs : [D] string;
  }
  
  // number of characters in common between two sequences
  proc ncommon(sequence1, sequence2) {
    var found = sequence1.length;
    if (sequence2.length < found) then
      found = sequence2.length;

    for (char1, char2, idx) in zip(sequence1(1..found), sequence2(1..found), 1..found) {
      if char1 != char2 {
        found = idx-1;
        break;
      }
    }
    return found;
  }

  // Naive suffix tree method: uses quadratic space, quadratic build time.
  // Inserts all suffixes in decreasing order of length
  class Edge {
    var parent: Node;
    var child: Node;
    var edgelabel: string;
    var len: int;

    proc splitAt(idx) {
      var newnode = new Node(parent.lab+edgelabel(1..idx), "");
      const newlab = edgelabel(idx+1..);
      var newedge = new Edge(newnode, child, newlab, newlab.length);
      newnode.add_out_edge(newedge);

      edgelabel = edgelabel(1..idx);
      len = edgelabel.length;
      child = newnode;

      return newnode;
    }
  }

  class Node {
    var edgekeys: domain(string);
    var edges: [edgekeys] Edge;
    var lab: string;
    var left_char: string;

    proc init(node_label: string, left: string) {
      lab = node_label;
      left_char = left;
    }

    proc edge(key) {
      return edges[key];
    }

    proc add_out_edge(edge) {
      const start = edge.edgelabel[1];
      edges[start] = edge;
    }

    // Walks the tree matching a string as far as it can.
    // Returns the last node found on the path, the first
    // character of the edge (if one exists), the index
    // into the edge, and the remaining string where the
    // match stops.
    proc walk(sequence) : (Node, string, int, string) {
      const start = sequence[1];

      if sequence == "" || !edgekeys.member(start) then
          return (this, "", 0, sequence);

      var edge = edges[start];
      const n = ncommon(sequence, edge.edgelabel);

      if n == edge.len then
        return edge.child.walk(sequence(n+1..));

      return (this, start, n, sequence(n+1..));
    }

    proc set_left_chars() {
      for edge in edges do
        edge.child.set_left_chars();

      const left_chars = [edge in edges] edge.child.left_char;
      if left_chars.size > 0 {
        left_char = left_chars[1];
        for i in 2..left_chars.size do 
          if left_chars[i] != left_char {
            left_char = "*";
            break;
          }
      }
    }

    proc printedges() {
      for key in edgekeys {
        writeln(edges[key].edgelabel);
        if edges[key].child != nil then
          edges[key].child.printedges();
      }
    }

    proc printnodes(prefix = "") {
      writeln(prefix, left_char, ' ', lab, ' ', edgekeys.size);
      for edge in edges do
        edge.child.printnodes(prefix + "  ");
    }

    proc maximal_repeats() : stringarray {
      var repeats : stringarray;
      const children = [edge in edges] edge.child;
      if children.size == 0 then
        return repeats;

      for child in children {
        var child_repeats = child.maximal_repeats();
        for child_repeat in child_repeats.strs do
          repeats.strs.push_back(child_repeat);
      }
      if left_char == "*" then
        repeats.strs.push_back(lab);

      return repeats;
    }
  }

  class QuadraticSuffixTree {
    var root: Node;

    proc insert(sequence, left="") {
      var (node, edgelab, ncommon, sremain) = root.walk(sequence);

      if sremain == "" then
        return;

      if edgelab != "" {
        var edge = node.edge(edgelab);
        node = edge.splitAt(ncommon);
      }

      var new_node = new Node(sequence, left);
      var newedge = new Edge(node, new_node, sremain, sremain.length);
      node.add_out_edge(newedge);
      return;
    }

    proc init(sequence) {
      root = new Node("", "");
      this.complete();
      var left = "";
      for i in 1..sequence.length {
        insert(sequence[i..], left);
        left = sequence[i];
      }
      root.set_left_chars();
    }

    proc printedges() {
      root.printedges();
    }

    proc printnodes() {
      root.printnodes();
    }

    proc maximal_repeats() {
      return root.maximal_repeats();
    }
  }

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

  proc main() {
    try! {
      for sequence in lines_from_file(infile, defaultfilename) {
        var suffixTree = new QuadraticSuffixTree(sequence+"$");
        const repeats = suffixTree.maximal_repeats();

        for rep in repeats.strs do
          if rep.length >= threshold then
            writeln(rep);

        delete suffixTree;
      }
    }
  }
}
