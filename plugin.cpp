#include <string.h>
#include <geoflow/gfSharedHeadersHash.h>
#include "register.hpp"

#define WIN_DECLSPEC
#ifdef WIN32
	#define WIN_DECLSPEC __declspec (dllexport)
#endif

extern "C"
{
	WIN_DECLSPEC geoflow::NodeRegister *allocator()
	{
    auto node_register = new geoflow::NodeRegister(GF_PLUGIN_NAME);
    register_nodes(*node_register);
		return node_register;
	}

	WIN_DECLSPEC void deleter(geoflow::NodeRegister *ptr)
	{
		delete ptr;
	}

	WIN_DECLSPEC void get_shared_headers_hash(char *hash)
	{
		strcpy(hash, GF_SHARED_HEADERS_HASH);
	}
}