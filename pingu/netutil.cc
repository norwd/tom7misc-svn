
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>
#include <netinet/in.h>
#include <arpa/inet.h>

#include <optional>
#include <string>
#include <cstdint>
#include <cstdio>
#include <cstring>

#include "netutil.h"

using namespace std;

std::string NetUtil::IPToString(uint32_t ip) {
  char buf[18] = {};
  sprintf(buf, "%d.%d.%d.%d",
		  (ip >> 24) & 255,
		  (ip >> 16) & 255,
		  (ip >> 8) & 255,
		  ip & 255);
  return buf;
}

std::optional<uint32_t> NetUtil::GetIPV4(const string &host,
										 string *canonical_name,
										 string *error) {
  struct addrinfo *ai_list;

  struct addrinfo ai_hints;
  memset (&ai_hints, 0, sizeof (ai_hints));
  ai_hints.ai_flags = 0;
  ai_hints.ai_flags |= AI_ADDRCONFIG;
  ai_hints.ai_flags |= AI_CANONNAME;
  ai_hints.ai_family = AF_INET;
  ai_hints.ai_socktype = SOCK_RAW;

  int ai_return = getaddrinfo (host.c_str(), nullptr, &ai_hints, &ai_list);
  if (ai_return != 0) {
	// Can get more info from return value
	if (error != nullptr) *error = "getaddrinfo failed";
	return {};
  }

  if (ai_list == nullptr) {
	if (error != nullptr) *error = "no hosts returned";
	return {};
  }

  // First entry gets canonical name.
  if (canonical_name != nullptr) {
	if (ai_list->ai_canonname != nullptr) {
	  *canonical_name = ai_list->ai_canonname;
	} else {
	  *canonical_name = host;
	}
  }

  for (struct addrinfo *ai_ptr = ai_list;
	   ai_ptr != nullptr;
	   ai_ptr = ai_ptr->ai_next) {
	// IPV4 only
	if (ai_ptr->ai_family == AF_INET) {
	  ai_ptr->ai_socktype = SOCK_RAW;
	  ai_ptr->ai_protocol = IPPROTO_ICMP;

	  #if 0
	  printf("AF_INET is %d\n"
			 "sock->sa_family is %d\n"
			 "addrlen is %d\n"
			 "sock->sa_data is ",
			 AF_INET,
			 ai_ptr->ai_addr->sa_family,
			 ai_ptr->ai_addrlen);
	  for (int i = 0; i < (int)ai_ptr->ai_addrlen; i++) {
		printf("%d.", (uint8_t)ai_ptr->ai_addr->sa_data[i]);
	  }

	  printf(" = ");
	  for (int i = 0; i < (int)ai_ptr->ai_addrlen; i++) {
		printf("%c", ai_ptr->ai_addr->sa_data[i]);
	  }
	  printf("\n");
	  #endif

	  sockaddr_in *s = (sockaddr_in*)ai_ptr->ai_addr;
	  const uint8_t *bytes = (const uint8_t*)&s->sin_addr.s_addr;
	  uint32_t ip = bytes[0] << 24 | bytes[1] << 16 | bytes[2] << 8 | bytes[3];
	  #if 0
	  printf("sin_family %d\n"
			 "sin_port %d\n"
			 "s_addr %08x\n"
			 "IP %s\n",
			 s->sin_family,
			 s->sin_port,
			 s->sin_addr.s_addr,
			 IPToString(ip).c_str());
	  #endif

	  freeaddrinfo(ai_list);
	  return {ip};
	}
  }

  if (error != nullptr) *error = "no AF_INET in list";
  freeaddrinfo(ai_list);
  return {};
}
