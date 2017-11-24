echo "[$(date +%T)] Running raytrace"
time ./bin/raytrace -r 720 -s 3 -o out/basic.png in/basic_pointlight/basic_pointlight.obj
time ./bin/raytrace -r 720 -s 3 -o out/simple.png in/simple_pointlight/simple_pointlight.obj
time ./bin/raytrace -r 720 -s 3 -o out/lines.png in/lines_pointlight/lines_pointlight.obj
time ./bin/raytrace -r 720 -s 3 -o out/instance.png in/instance10000_pointlight/instance10000_pointlight.obj
time ./bin/raytrace -r 720 -s 3 -o out/refl.png in/refl_pointlight/refl_pointlight.obj
time ./bin/raytrace -r 720 -s 3 -o out/transparentp.png in/transparentp_pointlight/transparentp_pointlight.obj
