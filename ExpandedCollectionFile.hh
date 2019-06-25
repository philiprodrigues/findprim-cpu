#ifndef EXPANDEDCOLLECTIONFILE_H
#define EXPANDEDCOLLECTIONFILE_H

#include "constants.hh"

#include <fstream>

class ExpandedCollectionFile
{
public:

    ExpandedCollectionFile(const char* filename)
        : m_file(filename, std::ifstream::binary),
          m_buffer(new char[collection_adcs_size])
    {
        if(m_file.bad() || m_file.fail() || !m_file.is_open()){
            throw std::runtime_error(std::string("Bad file ")+std::string(filename));
        }
        // Calculate the length of the file
        m_file.seekg(0, m_file.end);
        m_length = m_file.tellg();
        m_file.seekg(0, m_file.beg);
        if(m_length==0){
            throw std::runtime_error("Empty file");
        }
    }

    ~ExpandedCollectionFile()
    {
        m_file.close();
        delete[] m_buffer;
    }

    // Length of the file in bytes
    size_t length() const {return m_length;}
    // Number of messages in the file.
    size_t num_messages() const {return m_length/collection_adcs_size;}

    // Read the ith fragment into the buffer and return a pointer to
    // the first frame in the fragment. Subsequent calls will
    // overwrite the buffer with a different fragment
    MessageCollectionADCs* message(size_t i)
    {
        if(i>=num_messages()) return nullptr;
        // Seek to the right place in the file
        m_file.seekg(i*collection_adcs_size);
        // Check we didn't go past the end
        if(m_file.bad() || m_file.eof()) return nullptr;
        // Actually read the message into the buffer
        m_file.read(m_buffer,collection_adcs_size);
        if(m_file.bad() || m_file.eof()) return nullptr;
        return reinterpret_cast<MessageCollectionADCs*>(m_buffer);
    }

protected:
    std::ifstream m_file;
    size_t m_length;
    char* m_buffer;
};

/* Local Variables:  */
/* mode: c++         */
/* c-basic-offset: 4 */
/* End:              */

#endif
