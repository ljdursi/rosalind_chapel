module CSTR {
  use List;
  const defaultfilename="none";
  config const infile=defaultfilename;

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

  proc trivial(character) {
    const n = character.size;
    const s = + reduce(character);

    if (s == 0) || (s == 1) || (s == n-1) || (s == n) then
      return true;
    return false;
  }

  proc write_character(character) {
    for element in character do
      if element then
        write("1");
      else
        write("0");
    writeln();
  }

  proc character_table(lines) {
    const nlines = lines.size;
    const nnucleotides = lines[1].length;

    var table = new list([1..nlines] bool);
    for i in 1..nnucleotides {
      var character: [1..nlines] bool;
      character[1] = true;
      var onebase = lines[1][i];

      for j in 2..nlines do
        character[j] = (lines[j][i] == onebase);

      table.append(character);
    }

    return table;
  }
  
  proc main() {
    try! {
      var lines : [1..0] string;
      for line in lines_from_file(infile, defaultfilename) do
        lines.push_back(line);

      const table = character_table(lines);
      for character in table do
        if !trivial(character) then
          write_character(character);
    }
  }
}
