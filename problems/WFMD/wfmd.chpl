module WFMD {
  use LinearAlgebra;
  const defaultfilename="none";
  config const infile=defaultfilename;

  // use stdin if filename == defaultfilename
  proc params_from_file(filename: string, defaultfilename: string) throws {
    var n, m, g, k: int;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var success = channel.read(n, m, g, k);
    return (n, m, g, k);
  }

  proc n_choose_k(n: int, k: int): int {
    // pascal's triangle
    var triangle = Matrix(n+1, n+1, int);
    for i in 1..n+1 do
        triangle[i, 1] = 1;

    for i in 2..n+1 do
      for j in 2..i do
        triangle[i, j] = triangle[(i-1), (j-1)] + triangle[i-1, j];

    return triangle[n+1, k+1];
  }

  proc transition_matrix(npop) {
    const twon = 2*npop;
    var T = Matrix(twon+1, twon+1);  
    var twon_choose_k : [0..twon] real; 

    for i in 0..twon do
      twon_choose_k[i] = n_choose_k(twon, i);

    for j in 0..twon do
      for i in 0..twon do
        T[(j+1):int, (i+1):int] = twon_choose_k[j]*((1.0*(i)/twon)**j)*((1.0*(twon-(i))/twon)**(twon-(j)));

    return T;
  }

  proc evolve(initial_vec, npop, ngen) {
    var v = initial_vec;
    var T = transition_matrix(npop);

    for i in 1..ngen do
      v = dot(T, v);

    return v;
  }

  proc initial_vector(npop, m) {
    var v = Vector(2*npop+1);
    v[m+1] = 1;
    return v;
  }

  proc main() {
    try! {
      const (npop, m, ngen, k) = params_from_file(infile, defaultfilename);
      var v = initial_vector(npop, m);
      v = evolve(v, npop, ngen);

      writeln(+ reduce v[1..2*npop-k+1]);
    }
  }
}
