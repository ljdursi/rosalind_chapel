module FIB {
  const defaultfilename="none";
  config const infile=defaultfilename;

  // use stdin if filename == defaultfilename
  proc param_from_file(filename: string, defaultfilename: string) {
    var n: uint;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var success = channel.read(n);
    return n;
  }

  proc factorial(n:uint) : uint {
    if n == 0 then
      return 0;
    var result : uint = 1;
    for i in 2..n do
      result *= i;
    return result;
  }

  iter array_permutations(n: uint, arr : [1..n] uint): [1..n] uint {
    if n <= 1 then
      yield arr;
    else {
      var tail : [1..n-1] uint = arr(2..n);
      for perm in array_permutations( n-1, tail ) {
        var outperm = arr;
        for i in 1..n {
          outperm[1..i-1] = perm[1..i-1];
          outperm[i] = arr[1];
          outperm[i+1..] = perm[i..];
          yield outperm;
        }
      }
    }
  }

  proc main() {
    var n = param_from_file(infile, defaultfilename);
    var input : [1..n] uint = 1..n;
    writeln(factorial(n));

    for perm in array_permutations(n, input) do
        writeln(perm);
  }
}
