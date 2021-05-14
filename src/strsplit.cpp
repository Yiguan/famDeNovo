#include <boost/algorithm/string.hpp>
#include <vector>
void mysplit(std::vector<std::string> & sv, std::string & str, std::string sep)
{
	boost::split(sv, str, boost::is_any_of(sep));
}
