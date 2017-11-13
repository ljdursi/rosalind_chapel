#!/bin/bash
cd || exit
if [[ ! -f ~/chapel-1.16.0/util/setchplenv.bash ]]
then
    wget https://github.com/chapel-lang/chapel/releases/download/1.16.0/chapel-1.16.0-1.tar.gz
    tar xzf chapel-1.16.0-1.tar.gz
    rm chapel-1.16.0-1.tar.gz
    cd chapel-1.16.0 || exit
    source util/setchplenv.bash
    make
    export CHPL_AUX_FILESYS=curl
    patch -p0 <<EOF
--- modules/packages/Curl.chpl  2017-10-10 17:45:03.000000000 -0400
+++ modules/packages/Curl-new.chpl  2017-10-13 13:03:02.000000000 -0400
@@ -299,18 +299,18 @@
    :arg arg: the value to set the curl option specified by opt.
    :type arg: `int`, `string`, `bool`, or `slist`
 */
-proc file.setopt(opt:c_int, arg):bool {
+proc file.setopt(opt:c_int, arg):bool throws {
   var err:syserr = ENOERR;

   if (arg.type == slist) && (slist.home != this.home) {
-    ioerror(EFAULT:syserr, "in file.setopt(): slist, and curl handle do not reside on the same locale");
+    try ioerror(EFAULT:syserr, "in file.setopt(): slist, and curl handle do not reside on the same locale");
   }

   on this.home {
     err = chpl_curl_set_opt(this._file_internal, opt, arg);
   }

-  if err then ioerror(err, "in file.setopt(opt:c_int, arg)");
+  if err then try ioerror(err, "in file.setopt(opt:c_int, arg)");
   return true;
 }
EOF
    make
fi
