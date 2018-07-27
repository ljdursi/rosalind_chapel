module GREP {
  param defaultfilename="none";
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

  record arrnode {
    var D: domain(1);
    var l: [D] int;
  }
  record arrstring {
    var l: [1..0] string;
  }

  proc debruijn_graph(kmers) {
    const nkmers = kmers.size;
    var edges : [1..nkmers] domain(int);

    var prefixdomain: domain(string);
    var prefixes: [prefixdomain] domain(int);

    for (idx, kmer) in zip(1..nkmers, kmers) {
      const prefix = kmer[1..kmer.length-1];
      if !prefixdomain.member(prefix) then
          prefixdomain.add(prefix);
      prefixes[prefix].add(idx);
    }

    for (idx, kmer) in zip(1..nkmers, kmers) {
      const suffix = kmer(2..kmer.length);
      for neighbor in prefixes[suffix] do
        edges[idx].add(neighbor);
    }

    return edges;
  }

  proc dfs_paths(graph, start, sequences) {
    var pathstrings: domain(string);
    var stack: [1..0] (int, arrnode);

    var startpath : arrnode;
    startpath.l.push_back(start);

    stack.push_back((start, startpath));
    while stack.size > 0 {
        var (vtx, path) = stack[1];
        stack.remove(1);

        const neighbors = graph[vtx] - path.l : domain(int);
        if neighbors.size == 0 && path.l.size == sequences.size {
            var str = "";
            for idx in path.l do
                str += sequences[idx](1);

            pathstrings.add(str);
        }

        for nextnode in neighbors {
            var np = path;
            np.l.push_back(nextnode);
            stack.insert(1, (nextnode, np));
        }

    }
    return pathstrings;
  }

  proc main() {
    try! {
      var sequences : [1..0] string;
      for line in lines_from_file(infile, defaultfilename) do
        sequences.push_back(line);
     
      const nsequences = sequences.size;
      const graph = debruijn_graph(sequences);

      const paths = dfs_paths(graph, 1, sequences);
      for path in paths do
        writeln(path);
    }
  }
}
