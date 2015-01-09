#!/bin/bash
cur_dir=`dirname $0`
pkg_dir=${cur_dir}/../pkg/
bin_dir=${cur_dir}/../bin/
lib_dir=${cur_dir}/../lib/
mkdir -p ${bin_dir}
mkdir -p ${lib_dir}
echo "downloading Saifish..."
salmon_url=https://github.com/kingsfordgroup/sailfish/releases/download/v0.2.3/SalmonBeta-v0.2.3_ubuntu-14.04.tar.gz
salmon_tarball=${pkg_dir}${salmon_url##*/}
wget -P ${pkg_dir} -N https://github.com/kingsfordgroup/sailfish/releases/download/v0.2.3/SalmonBeta-v0.2.3_ubuntu-14.04.tar.gz
tar -zxvf ${salmon_tarball} -C ${pkg_dir}
echo "copying Sailfish executables and libraries to bin/ lib/"
salmon_dir=${salmon_tarball%.tar.gz}
cp ${salmon_dir}/bin/* ${bin_dir}
cp ${salmon_dir}/lib/* ${lib_dir}
echo "downloading Star..."
star_url=https://github.com/alexdobin/STAR/archive/STAR_2.4.0h1.tar.gz
star_tarball=${pkg_dir}${star_url##*/}
wget -P ${pkg_dir} -N ${star_url}
tar -zxvf ${star_tarball} -C ${pkg_dir}
echo "copying star executable to bin/"
star_dir=${pkg_dir}STAR-STAR_2.4.0h1
cp ${star_dir}/bin/Linux_x86_64/* ${bin_dir}
echo "done"

