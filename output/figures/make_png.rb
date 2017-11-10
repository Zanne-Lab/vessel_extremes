Dir["*.pdf"].each do |file|
  if file.include? "pdf"
    puts file
    `convert -density 600 #{file} #{file.sub("pdf","png")}`
  end
end
