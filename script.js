document.addEventListener("DOMContentLoaded", function() {
  const numberOfRaindrops = 1; // Adjust the number of raindrops
  const raindrops = document.getElementById("raindrops");

  for (let i = 0; i < numberOfRaindrops; i++) {
    const raindrop = document.createElement("div");
    raindrop.className = "raindrop";
    raindrops.appendChild(raindrop);

  }
});
