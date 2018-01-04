module PRSM {
  use Sort;
  use Search;
  use List;
  const defaultfilename = "none";
  config const infile = defaultfilename;
  config const weightfile = "weights.txt";
  config const tol = 5.e-5;

  proc readtable(weightfile) throws {
    var aminos: [1..0] string;
    var weights: [1..0] real;

    try! {
      var file = open(weightfile, iomode.r);
      var channel = file.reader();

      var amino : string;
      var weight : real;
      while (channel.readln(amino, weight)) {
        aminos.push_back(amino);
        weights.push_back(weight);
      }
    }

    return (aminos, weights);
  }

  proc reals_from_file(filename: string, defaultfilename: string) throws {
    var result : [1..0] real;
    var line: string;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    while (channel.readline(line)) do 
      result.push_back(line: real);

    return result;
  }

  proc nearest(sortednums, query) {
    var (found, idx) = binarySearch(sortednums, query);
    var delta = 100000.0;
    var bestidx = -1;

    for posidx in max(1, idx-1)..min(idx+1, sortednums.size) {
      var thisdelta = abs(sortednums[posidx]-query);
      if delta < thisdelta then
        continue;

      delta = thisdelta;
      bestidx = posidx;
    }

    return (bestidx, delta);
  }

  proc spectrum_graph(peptides, aminos, weights, tol=5.e-5) {
    const np = peptides.size;
    const na = aminos.size;
    var graph : [1..np] list((string, int));

    for i in 1..np {
      for j in 1..na {
        const neighweight = peptides[i] + weights[j];
        const (neigh, diff) = nearest(peptides, neighweight);
        if diff < tol then
          graph[i].append((aminos[j], neigh));
      }
    }

    return graph;
  }

  proc longest_path(graph, peptides) {
    const n = peptides.size;
    var best_pathlen : [1..n] int;
    var last : [1..n] int;
    var aminos : [1..n] string;

    aminos = '-';
    last = -1;
    best_pathlen = 0;

    var longest = 0;
    var longest_end = -1;

    for (i, peptide) in zip(1..n, peptides) {
      const curbest = best_pathlen[i];
      for (amino, neighbor) in graph[i] {
        const current = curbest + 1;
        if current > best_pathlen[neighbor] {
          best_pathlen[neighbor] = current;
          last[neighbor] = i;
          aminos[neighbor] = amino;
          if (current > longest) {
            longest = current;
            longest_end = neighbor;
          }
        }
      }
    }

    var longest_protein = "";
    var idx = longest_end;
    while (last[idx] > -1) {
      longest_protein = aminos[idx] + longest_protein;
      idx = last[idx];
    }

    return longest_protein;
  }

  proc main() {
    try! {
      var (aminos, weights) = readtable(weightfile);
      var peptides = reals_from_file(infile, defaultfilename);
      sort(peptides);
      var graph = spectrum_graph(peptides, aminos, weights);
      var result = longest_path(graph, peptides);
      writeln(result);
    }
  }
}
