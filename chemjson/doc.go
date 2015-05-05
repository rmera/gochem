package chemjson

//Package chemjson implements serializacion and unserialization of
//goChem data types. It's planned use is the communication of goChem
//programs with other, independent programs which can be written in
//languages other than Go, as long as those languages implement a
//way of serializing and unserializing JSON data and some library
//to deal with the goChem types is implemented.
//chemjson also implements the transmision of options, so an external
//program can transmit data an options for a job to a goChem program
//and later collect the results, for instance, via UNIX pipes.
