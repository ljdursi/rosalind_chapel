module QRT {
  const defaultfilename="none";
  config const infile=defaultfilename;

  // use stdin if filename == defaultfilename
  proc lines_from_file(filename: string, defaultfilename: string) throws {
    var channel = stdin;

    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var lines: [1..0] string;
    var line: string;
    while (channel.readline(line)) do
      lines.push_back(line.strip());

    return lines;
  }

  proc partition_from_character(character: string) {
    var part_one : [1..0] int;
    var part_two : [1..0] int;

    for (i, c) in zip(1.., character) do
      select(c) {
        when("0") do
          part_one.push_back(i);
        when("1") do
          part_two.push_back(i);
      } 

    return (part_one, part_two);
  }

  iter pairs(inputlist) {
    const lo = inputlist.domain.dim(1).low;
    const hi = inputlist.domain.dim(1).high;

    for i in lo..hi do
      for j in i+1..hi do
        yield (inputlist[i], inputlist[j]);
  }

  iter quartets_from_partition(part_one, part_two) {
    for pair1 in pairs(part_one) do
      for pair2 in pairs(part_two) do
        if pair1[1] <= pair2[1] then
          yield (pair1[1], pair1[2], pair2[1], pair2[2]);
        else
          yield (pair2[1], pair2[2], pair1[1], pair1[2]);
  }

  proc quartet_string(allnames, quartet) {
    var result = "{" + allnames[quartet[1]] + ", " + allnames[quartet[2]] + "}";
    result += " {" + allnames[quartet[3]] + ", " + allnames[quartet[4]] + "}";

    return result;
  }

  proc main() {
    try! {
      var lines = lines_from_file(infile, defaultfilename);
      const allnames = lines[1].split();

      var quartets : domain(4*int);
      for line in lines[2..] {
        var (part1, part2) = partition_from_character(line);
        for quartet in quartets_from_partition(part1, part2) do
          quartets.add(quartet);
      }

      for quartet in quartets do
        writeln(quartet_string(allnames, quartet));
    }
  }
}
