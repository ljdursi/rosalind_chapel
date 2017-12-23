module PRSM {
  use Sort;
  const defaultfilename = "none";
  config const infile = defaultfilename;
  config const weightfile = "weights.txt";
  config const tol = 1.e-7;

  proc readtable(weightfile) throws {
    var aminonames: domain(string);
    var amino_weights : [aminonames] real;

    try! {
      var file = open(weightfile, iomode.r);
      var channel = file.reader();

      var amino : string;
      var weight : real;
      while (channel.readln(amino, weight)) {
        aminonames.add(amino);
        amino_weights[amino] = weight;
      }
    }

    return amino_weights;
  }

  proc lines_from_file(filename: string, defaultfilename: string) throws {
    var result : [1..0] string;
    var line: string;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    while (channel.readline(line)) do 
      result.push_back(line.strip());

    return result;
  }

  proc spectrum(sequence, aminoweights) {
    const n = sequence.length;
    var weights: [1..n] real;
    for (i, amino) in zip(1..n, sequence) do
      weights[i] = aminoweights[amino];

    const totalweight = +reduce(weights);
    const prefixes = +scan(weights);
    const suffixes = totalweight - prefixes[1..n-1];

    var spec : [1..2*n] real;
    spec[1..n] = prefixes;
    spec[n+1..2*n-1] = suffixes;
    spec[2*n] = totalweight;
    sort(spec);

    return spec;
  }

  proc differences(a, b) {
    var diffs : [1..0] real;
    for ia in a do
      for ib in b do
        diffs.push_back(ia-ib);
    
    sort(diffs);
    return diffs;
  }

  proc max_multiplicity(sortednums, tol) {
    var current_mult = 0;
    var current_last = -1.0;
    var best_mult = 0;

    for item in sortednums {
      if abs(item - current_last) < tol then
        current_mult += 1;
      else {
        if current_mult > best_mult then
          best_mult = current_mult;
        current_last = item;
        current_mult = 1;
      }
    }

    return max(current_mult, best_mult);
  }

  proc main() {
    try! {
      var amino_weights = readtable(weightfile);
      var lines = lines_from_file(infile, defaultfilename);
      const m = lines.size;

      const n = lines[1]:int;
      var peaks : [1..m-n-1] real;
      for (i, line) in zip(1.., lines[n+2..]) do
        peaks[i] = line: real;
      sort(peaks);

      var best = (-1, "");
      for protein in lines[2..n+1] {
        const spec = spectrum(protein, amino_weights);
        const diffs = differences(spec, peaks);
        const maxmult = max_multiplicity(diffs, tol);
        if (maxmult, protein) > best then
          best = (maxmult, protein);
      }

      writeln(best[1]);
      writeln(best[2]);
    }
  }
}
