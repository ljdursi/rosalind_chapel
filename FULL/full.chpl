module FULL {
  use Search;
  use Sort;
  const defaultfilename = "none";
  config const infile = defaultfilename;
  config const weightfile = "weights.txt";
  config const tol = 1.e-7;

  proc readtable() throws {
    var rows : [1..0] (real, string);

    try! {
      var file = open(weightfile, iomode.r);
      var channel = file.reader();

      var amino : string;
      var weight : real;
      while (channel.readln(amino, weight)) do
        rows.push_back((weight, amino));
    }
    sort(rows);

    var weights : [1..rows.size] real;
    var aminos: [1..rows.size] string;
    var idx = 1;
    for (weight, amino) in rows {
        weights[idx] = weight;
        aminos[idx] = amino;
        idx += 1;
    }
    return (aminos, weights);
  }

  // use stdin if filename == defaultfilename
  proc reals_from_file(filename: string, defaultfilename: string) throws {
    var result : [1..0] real;
    var item: real;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    while (channel.read(item)) do 
      result.push_back(item);

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

  proc infer_protein(prefixes, aminos, weights, tol) {
    var protein = "";

    var (startweight, endweight) = (prefixes[1], prefixes[prefixes.size]);
    prefixes.remove(prefixes.size);
    prefixes.remove(1);
    const n = prefixes.size / 2;

    for loc in 1..n {
      var best = (1000.0, -1, -1, -1, 0.0);

      for (j, weight) in zip(1.., weights) {
        const (newstartweight, newendweight) = (startweight+weight, endweight-weight);
        const (test_prefix, delta_prefix) = nearest(prefixes, newstartweight);
        const (test_suffix, delta_suffix) = nearest(prefixes, newendweight);

        const delta = max(delta_prefix, delta_suffix);
        if delta < best[1] then
          best = (delta, j, test_prefix, test_suffix, weight);

        if delta < tol then
          break;
      }

      const (delta, j, test_prefix, test_suffix, w) = best;
      protein = protein + aminos[j];
      startweight += w;
      endweight -= w;

      prefixes.remove(max(test_prefix, test_suffix));
      prefixes.remove(min(test_prefix, test_suffix));
    }
    return protein;
  }

  proc main() {
    try! {
      var (aminos, weights) = readtable();
      var prefixes = reals_from_file(infile, defaultfilename);
      prefixes.remove(1);

      writeln(infer_protein(prefixes, aminos, weights, tol));
    }
  }
}
