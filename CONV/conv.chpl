module CONV {
  use Sort;
  const defaultfilename="none";
  config const infile=defaultfilename;

  proc strings_from_file(filename: string, defaultfilename: string) throws {
    var text: string = "";
    var line: string;
    var lines: [1..0] string;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    while (channel.readline(line)) {
      lines.push_back(line.strip());
    }

    return lines;
  }

  proc differences(a, b, delta=1.0e-6) {
    var diffs :[1..(a.size)*(b.size)] real;

    var i = 1;
    for aa in a do
      for bb in b {
        diffs[i] = aa - bb;
        i += 1;
      }

    var multiplicities: [1..0] (uint, real);
    var last = -999999999.0;
    var count = 0 : uint;
    for item in sorted(diffs) {
      if (abs(item-last) < delta) then
        count += 1;
      else {
        multiplicities.push_back((count, last));
        last = item;
        count = 1;
      }
    }
    multiplicities.push_back((count, last));
    multiplicities.pop_front();

    return multiplicities;
  }

  proc main() {
    try! {
      var lines = strings_from_file(infile, defaultfilename);
      var a : [1..0] real;
      var b : [1..0] real;

      for item in lines[1].split() do
        a.push_back(item:real);
      for item in lines[2].split() do
        b.push_back(item:real);

      var diffs = differences(a, b);
      sort(diffs);

      const lastdiff = diffs[diffs.size+1];
      writeln(lastdiff[1]);
      writeln(abs(lastdiff[2]));
    }
  }
}
