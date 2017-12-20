module PCOV {
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

  record arrstring {
    var l: [1..0] string;
  }

  proc debruijn_graph(kmers) {
    var prefixdomain: domain(string);
    var prefixes: [prefixdomain] arrstring;
    var edgedomain: domain(string);

    for kmer in kmers {
      edgedomain.add(kmer);
      const prefix = kmer[1..kmer.length-1];
      prefixdomain.add(prefix);
    }
    for kmer in kmers{
      const prefix = kmer[1..kmer.length-1];
      prefixes[prefix].l.push_back(kmer);
    }

    var edges : [edgedomain] arrstring;
    for kmer in kmers {
      const suffix = kmer(2..kmer.length);
      for neighbor in prefixes[suffix].l do
        edges[kmer].l.push_back(neighbor);
    }

    return edges;
  }

  proc reconstruct_chromosome(edges) {
    const nodes = edges.domain.sorted();
    var current = nodes[1];
    var chromosome = current;
    const k = nodes[1].length;

    for i in nodes {
       current = edges[current].l[1];
       chromosome += current(k);
    }

    return chromosome(1..chromosome.length-k);
  }

  proc main() {
    try! {
      var sequences : [1..0] string;
      for line in lines_from_file(infile, defaultfilename) do
        sequences.push_back(line);
     
      var graph = debruijn_graph(sequences);
      const chromosome = reconstruct_chromosome(graph);
      writeln(chromosome);
    }
  }
}
