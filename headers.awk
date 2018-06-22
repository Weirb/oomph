START
{
 print "// This file was generated automatically during the make process"
 print "// and it will be remade automatically"
}
{
 # Loop over all header file names
 for (i=1;i<=NF;i++) printf "#include<"$i"> \n"
} 


