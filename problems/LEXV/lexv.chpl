module DNA {
  const defaultfilename="none";
  config const infile=defaultfilename;

  // use stdin if filename == defaultfilename
  proc params_from_file(filename: string, defaultfilename: string) throws {
    var alphabetstring: string = "";
    var k: uint;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var success = channel.readline(alphabetstring);
    success = channel.read(k);

    var alphabet: [1..0] string;
    for base in alphabetstring.split() do
        alphabet.push_back(base);

    return (alphabet, k);
  }

  iter lex_strings_from_alphabet(k:uint, n:int, alphabet:[1..n] string, first: bool): string {
    if k == 0 then
        yield "";
    else {
        if !first then
            yield "";
        for base in alphabet do
            for tail in lex_strings_from_alphabet(k-1, n, alphabet, false) do
                yield base + tail;
    }
  }

  proc main() {
    try! {
      var params = params_from_file(infile, defaultfilename);
      var alphabet = params[1];
      var k : uint = params[2];
      var n : int = alphabet.size;

      for kmer in lex_strings_from_alphabet(k, n, alphabet, true) do
          writeln(kmer);
    }
  }
}
