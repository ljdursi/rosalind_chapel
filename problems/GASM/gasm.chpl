module GASM {
  // Genome Assembly With Reads
  // http://rosalind.info/problems/gasm/
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

  proc reverse_complement(seq: string) {
    // Reverse complement of a string
    const complement = ['A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A'];

    var rc = "";
    for i in 1..seq.length by -1 do
      rc = rc + complement[seq(i)];

    return rc;
  }

  proc cyclic_rc_match(seq1: string, seq2: string) {
    // Test to see if two cyclic strings are consistent with
    // being a RC of each other
    if seq1.length != seq2.length then
      return false;

    const rc_seq2 = reverse_complement(seq2);
    const doubled_seq1 = seq1 + seq1;

    const loc = doubled_seq1.find(rc_seq2);
    return loc != 0;
  }

  proc generate_kmers(reads, k) {
    // List of k-mers of consistent with the
    // reads, and their reverse complements
    var kmers : [1..0] string;
    const n = reads[1].length;
    for read in reads {
      for i in 1..n-k+1 do
        kmers.push_back(read(i..i+k-1));
      const rc_read = reverse_complement(read);
      for i in 1..n-k+1 do
        kmers.push_back(rc_read(i..i+k-1));
    }
    return kmers;
  }

  // we need this type definition so that we can 
  // have assiciative arrays of them
  // Note that if you don't include the domain (and make it first!),
  // you'll need an explicit copy constructor to make operations on
  // arrays of these things not fail
  record arrstring {
    var D: domain(1);
    var l : [D] string;
  }

  proc debruijn_graph(kmers) {
    // build a deBruijn graph of the kmers
    // returns associative array: string -> arrstring
    //                            in_vtx -> [outvtx]
    var prefixdomain: domain(string);
    var prefixes: [prefixdomain] domain(string);
    var edgedomain: domain(string);
    var edges : [edgedomain] domain(string);

    for kmer in kmers {
      edgedomain.add(kmer);
      const prefix = kmer[1..kmer.length-1];
      prefixdomain.add(prefix);
      prefixes[prefix].add(kmer);
    }

    for kmer in kmers {
      const suffix = kmer(2..kmer.length);
      if !prefixdomain.member(suffix) then
        continue;
      for neighbor in prefixes[suffix] do
        edges[kmer].add(neighbor);
    }

    return edges;
  }

  proc get_next(ref graph, key) {
    if !graph.domain.member(key) then
      return "";

    var next = "";
    if !graph[key].isEmpty() then
      for item in graph[key] {
        next = item;
        graph[key].remove(item);
        break;
      }
    
    if graph[key].isEmpty() then
      graph.domain.remove(key);
    
    return next;
  }

  proc reconstruct_chromosome(edges) {
    var chromosomes : [1..0] string;

    while edges.domain.size > 0 {
      const nodes = edges.domain.sorted();
      var current = nodes[1];
      const k = current.length;

      var cur_chromosome = current;
      var nextnode = get_next(edges, current);
      while nextnode != "" {
        cur_chromosome = cur_chromosome + nextnode(k);
        current = nextnode;
        nextnode = get_next(edges, nextnode);
      }

      const totlen = cur_chromosome.length;
      chromosomes.push_back(cur_chromosome(1..totlen-k));
    }

    return chromosomes;
  }

  proc success(chromosomes) {
    if chromosomes.size != 2 then
      return false;
    return cyclic_rc_match(chromosomes[1], chromosomes[2]);
  }

  proc main() {
    try! {
      var sequences : [1..0] string;
      for line in lines_from_file(infile, defaultfilename) do
        sequences.push_back(line);
     
      for k in 2..sequences[1].length by -1 {
        const kmers = generate_kmers(sequences, k);
        var graph = debruijn_graph(kmers);
        const chromosomes = reconstruct_chromosome(graph);
        if success(chromosomes) {
          writeln(chromosomes[1]);
          break;
        }
      }
    }
  }
}
