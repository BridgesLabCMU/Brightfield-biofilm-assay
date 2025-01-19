module utils

export unzip_folder, write_zip

using ZipFile

function unzip_folder(zip_path)
    # Destination folder name derived from zip path
    destination_path = replace(zip_path, r"\.zip$" => "")

    # Create the destination folder if it doesn't exist
    if !isdir(destination_path)
        mkdir(destination_path)
    end

    # Open the zip file
    r = ZipFile.Reader(zip_path)
    try
        for file in r.files
            # Construct the full path for the destination file
            file_destination = joinpath(destination_path, file.name)

            # Check if the entry is a directory
            if endswith(file.name, "/")
                if !isdir(file_destination)
                    mkdir(file_destination)
                end
            else
                # Ensure parent directories exist
                parent_dir = dirname(file_destination)
                if !isdir(parent_dir)
                    mkdir(parent_dir)
                end

                # Extract the file by directly reading its contents
                open(file_destination, "w") do output_stream
                    write(output_stream, read(file))
                end
            end
        end
    finally
        close(r)  # Ensure the ZipFile.Reader is closed
    end
end

function write_zip(DOWNLOADS_ZIP, DOWNLOADS_PATH)
	compress = true
	zdir = ZipFile.Writer(DOWNLOADS_ZIP) 
	for (root, dirs, files) in walkdir(DOWNLOADS_PATH)
		for file in files
			filepath = joinpath(root, file)
			f = open(filepath, "r")
			content = read(f, String)
			close(f)

			zippath = relpath(filepath, DOWNLOADS_PATH * "/")
			println("\t$(zippath)")
			zf = ZipFile.addfile(zdir, zippath; method=(compress ? ZipFile.Deflate : DOWNLOADS_ZIP.Store));
			write(zf, content)
		end
	end
	close(zdir)
end

end
