module NWCK {
  use Regexp;
  const defaultfilename="none";
  config const infile=defaultfilename;

  proc params_from_file(filename: string, defaultfilename: string) throws {
    var line: string;
    var n : uint;
    var line1, line2: string;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var success = channel.readln(n);
    success = channel.readline(line1);
    success = channel.readline(line2);

    return (n, line1, line2);
  }

  proc set_from_string(line) : domain(uint) throws {
    var result : domain(uint);
    var splitter = compile("{|}|,|[[:space:]]+");
    for item in splitter.split(line) {
      if item != "" {
        result.add(item:uint);
      }
    }

    return result;
  }


  proc main() {
    try! {
      var (n, aline, bline) = params_from_file(infile, defaultfilename);
      var a = set_from_string(aline);
      var b = set_from_string(bline);
      var u : domain(uint) = 1..n;

      writeln(a + b);
      writeln(a & b);
      writeln(a - b);
      writeln(b - a);
      writeln(u - a);
      writeln(u - b);
    }
  }
}
