/**
 * \file FilesystemUtils.h
 *
 * \brief Filesystem utilities.
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 4/7/11
 */

#ifndef FILESYSTEMUTILS_H
#define	FILESYSTEMUTILS_H

#include <string>

/***************************************************************************//**
 * Get the full filename without the extension.
 * \param [in] fullFilename complete filename
 * \return path/filename without extension
 ******************************************************************************/
std::string GetFileBasename(std::string fullFilename);

/***************************************************************************//**
 * Get the filename extension.
 * \param [in] fullFilename complete filename
 * \return filename extension
 ******************************************************************************/
std::string GetFileExtension(std::string fullFilename);

#endif	/* FILESYSTEMUTILS_H */

