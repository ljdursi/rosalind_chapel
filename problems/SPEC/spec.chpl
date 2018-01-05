module RNA {
  use Search;
  use Sort;
  const defaultfilename="none";
  config const infile=defaultfilename;
  config const weightfile="weights.txt";

  proc readtable() throws {
    var rows : [1..0] (real, string);
    rows.push_back((-1.0, "None"));

    try! {
      var file = open(weightfile, iomode.r);
      var channel = file.reader();

      var amino : string;
      var weight : real;
      while (channel.readln(amino, weight)) do
        rows.push_back((weight, amino));
    }
    rows.push_back((99999.0, "None"));
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

    while (channel.readln(item)) do 
      result.push_back(item);

    return result;
  }

  proc infer_protein(prefixes, aminos, weights) {
    var protein = "";
    sort(prefixes);

    for idx in 2..prefixes.size {
      const delta = prefixes[idx] - prefixes[idx-1];
      var (found, loc) = binarySearch(weights, delta);
      if abs(weights[loc-1]-delta) < abs(weights[loc]-delta) then
        loc = loc-1;
      protein = protein + aminos[loc];
    }
    return protein;
  }

  proc main() {
    try! {
      var (aminos, weights) = readtable();
      var prefixes = reals_from_file(infile, defaultfilename);
      
      writeln(infer_protein(prefixes, aminos, weights));
    }
  }
}
