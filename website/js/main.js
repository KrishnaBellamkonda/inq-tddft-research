// Minimal progressive enhancement for the static project site.
// 1. Smooth-scroll is handled in CSS (html { scroll-behavior: smooth }).
// 2. Reveal sections on scroll for a little polish.
// 3. Pause the hero video when off-screen (saves battery, keeps it tidy).

document.addEventListener("DOMContentLoaded", () => {
  // Fade sections in as they enter the viewport.
  const blocks = document.querySelectorAll("section.block, details.appendix");
  if ("IntersectionObserver" in window) {
    blocks.forEach((b) => {
      b.style.opacity = "0";
      b.style.transform = "translateY(16px)";
      b.style.transition = "opacity 0.6s ease, transform 0.6s ease";
    });
    const io = new IntersectionObserver(
      (entries) => {
        entries.forEach((e) => {
          if (e.isIntersecting) {
            e.target.style.opacity = "1";
            e.target.style.transform = "none";
            io.unobserve(e.target);
          }
        });
      },
      { threshold: 0.12 }
    );
    blocks.forEach((b) => io.observe(b));
  }

  // Pause hero video when scrolled past (it's decorative once you're reading).
  const video = document.querySelector(".hero video");
  if (video && "IntersectionObserver" in window) {
    const vio = new IntersectionObserver(
      (entries) => {
        entries.forEach((e) => {
          if (e.isIntersecting) video.play().catch(() => {});
          else video.pause();
        });
      },
      { threshold: 0.25 }
    );
    vio.observe(video);
  }
});
