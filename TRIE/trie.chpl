module GRPH {
  const defaultfilename="none";
  config const infile=defaultfilename;

  class Node {
    var num: int;
    var bases: domain(string);
    var edges: [bases] Node;
  }

  class Trie {
    var nnodes : int = 1;
    var root = new Node(1);

    proc insert(sequence: string) {
      var node = root;
      for base in sequence {
        if node.bases.member(base) then
          node = node.edges[base];
        else {
          nnodes += 1;
          node.edges[base] = new Node(nnodes);
          node = node.edges[base];
        }
      }
    }

    proc edgelist() {
      proc dfs_edgelist(n: Node) {
        for base in n.bases {
          writeln(n.num, " ", n.edges[base].num, " ", base);
          dfs_edgelist(n.edges[base]);
        }
      }
      dfs_edgelist(root);
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
    var trie = new Trie();
    try! {
      for sequence in lines_from_file(infile, defaultfilename) do
        trie.insert(sequence);
      }
    trie.edgelist();
  }
}
