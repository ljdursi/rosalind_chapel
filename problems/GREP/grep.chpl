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

  proc string_from_path(path, seqs) {
    var str = "";
    for idx in path.l do
        str += seqs[idx](1);
    
    return str;
  }

  proc dfs_paths(graph, start, sequences) {
    var pathstrings: domain(string);
    var stack: [1..0] arrnode;

    var startpath : arrnode;
    startpath.D = 1..1;
    startpath.l = [ start ];
    stack.push_back(startpath);

    while stack.size > 0 {
      const path = stack[stack.size];
      const vtx = path.l[path.l.size];
      stack.pop_back();

      const neighbors = graph[vtx] - path.l : domain(int);
      for nextnode in neighbors {
        var np = path;
        np.l.push_back(nextnode);
        stack.push_back(np);
      }

      if path.l.size == sequences.size then
        pathstrings.add( string_from_path(path, sequences) );

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
