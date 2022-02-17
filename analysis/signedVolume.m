function v = signedVolume(a, b, c)
v = sum(a.*cross(b, c));
end