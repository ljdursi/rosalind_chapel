module SIGN {
  const defaultfilename="none";
  config const infile=defaultfilename;

  // use stdin if filename == defaultfilename
  proc param_from_file(filename: string, defaultfilename: string) throws {
    var n: int;

    try! {
      var channel = stdin;
      if infile != defaultfilename {
        var file = open(infile, iomode.r);
        channel = file.reader();
      }
      var success = channel.read(n);
    }
    return n;
  }

  proc factorial(n:int) : int {
    if n == 0 then
      return 0;
    var result : int = 1;
    for i in 2..n do
      result *= i;
    return result;
  }

  iter array_permutations(n: int, arr : [1..n] int): [1..n] int {
    if n <= 1 then
      yield arr;
    else {
      var tail : [1..n-1] int = arr(2..n);
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

  iter signvecs(n: int): [1..n] int {
    var signs : [1..n] int = 1;
    for signval in 1..2**n {
      var result = signs;
      var s = signval;
      for j in 1..n {
        if s % 2 == 1 then
          result[j] = -1;
        s /= 2;
      }
      yield result;
    }
  }

  proc main() {
    try! {
      const n = param_from_file(infile, defaultfilename);
      const nsigns = 2**n;
      writeln(factorial(n) * nsigns);
  
      var input : [1..n] int = 1..n;
      for perm in array_permutations(n, input) {
        for signvec in signvecs(n) {
          for (item, sign) in zip(perm, signvec) do
            write(item*sign, " ");
          writeln();
        }
      }
    }
  }
}
