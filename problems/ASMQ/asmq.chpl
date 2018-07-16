module ASMQ {
  // Assessing Assembly Quality with N50 and N75
  // http://rosalind.info/problems/asmq/
  use Sort;
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

  proc main() {
    try! {
      var read_lengths : [1..0] uint;
      for line in lines_from_file(infile, defaultfilename) do
        read_lengths.push_back(line.length: uint);
     
      const n_tot : uint = + reduce read_lengths;
      sort(read_lengths);
      var tot_so_far : uint = 0;
      var n50, n75: uint;
      for rlen in read_lengths {
        if (n_tot - tot_so_far) >= 3*n_tot/4 then
          n75 = rlen;
        if (n_tot - tot_so_far) >= n_tot/2 then 
          n50 = rlen;
        tot_so_far += rlen;
      }

      writeln(n50, ' ', n75);
    }
  }
}
