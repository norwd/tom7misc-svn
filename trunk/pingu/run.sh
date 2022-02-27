#!/bin/bash

if [[ $(id -u) -ne 0 ]] ;
then echo "Please run as root" ;
	 exit -1 ;
fi

set +x

# TODO make this a command-line option
PLUGIN=compu
# PLUGIN=tetru
# PLUGIN=pingu
MOUNTPOINT="/mnt/$PLUGIN"
SOCKET="/tmp/nbdsocket.$PLUGIN"

# make sure this is loaded (e.g. after reboot)
modprobe nbd

# clean up any existing
umount "$MOUNTPOINT"
mkdir -p "$MOUNTPOINT"
nbd-client -d /dev/nbd0
killall -9 "$PLUGIN-viz.exe"
killall -9 nbdkit
rm -f "$SOCKET"

if [ "$1" = "stop" ]; then
	echo "Stopped."
	df -k "$MOUNTPOINT"
	exit 0
fi

# Note: 65536 bytes (128 blocks) is too small for FAT
# 

# note that the argument to --run includes an escaped $unixsocket;
# this is a nbdkit concept, not a bash variable.

../../nbdkit/server/nbdkit --verbose -U "$SOCKET" "./$PLUGIN.so" 51200 --run "./mount.sh \$unixsocket $MOUNTPOINT"
# 2>&1 | "viz/$PLUGIN-viz.exe"

# grep -v 'TVIZ\[r '
# drop-in replacement with memory plugin:
# ../../nbdkit/server/nbdkit --verbose -U "$SOCKET" ../../nbdkit/plugins/memory/.libs/nbdkit-memory-plugin.so 51200 --run "./mount.sh \$unixsocket $MOUNTPOINT"

