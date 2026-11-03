package com.nicksard.appsleep;

import android.accessibilityservice.AccessibilityService;
import android.content.Intent;
import android.content.pm.PackageManager;
import android.content.pm.ResolveInfo;
import android.os.Handler;
import android.os.Looper;
import android.os.SystemClock;
import android.view.accessibility.AccessibilityEvent;
import android.view.inputmethod.InputMethodInfo;
import android.view.inputmethod.InputMethodManager;

import java.time.LocalDateTime;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Watches which app comes to the foreground. If that app is asleep on an active schedule,
 * it sends the user to the home screen and shows a short "this app is asleep" screen.
 *
 * Android does not let ordinary apps freeze or force-stop other apps, so this is the same
 * approach app blockers on the Play Store use.
 */
public class SleepGuardService extends AccessibilityService {
    private static final long RECHECK_MS = 30_000;

    private final Handler handler = new Handler(Looper.getMainLooper());
    private final Set<String> passThrough = new HashSet<>();
    private String foregroundPackage;
    private String lastBlockedPackage;
    private long lastBlockedAt;

    // Catches the case where a schedule starts while you are already inside the app.
    private final Runnable recheck = new Runnable() {
        @Override
        public void run() {
            if (foregroundPackage != null) check(foregroundPackage);
            handler.postDelayed(this, RECHECK_MS);
        }
    };

    @Override
    protected void onServiceConnected() {
        super.onServiceConnected();
        refreshPassThrough();
        handler.removeCallbacks(recheck);
        handler.postDelayed(recheck, RECHECK_MS);
    }

    @Override
    public void onAccessibilityEvent(AccessibilityEvent event) {
        if (event == null || event.getEventType() != AccessibilityEvent.TYPE_WINDOW_STATE_CHANGED) {
            return;
        }
        CharSequence pkg = event.getPackageName();
        if (pkg == null) return;
        String p = pkg.toString();
        // The keyboard and system UI (notification shade, volume panel) pop up on top of
        // other apps; they don't tell us the foreground app changed.
        if (p.equals("com.android.systemui") || isInputMethod(p)) return;
        foregroundPackage = p;
        if (passThrough.contains(p)) return;
        check(p);
    }

    private void check(String p) {
        if (passThrough.contains(p)) return;
        Schedule s = ScheduleStore.sleepingScheduleFor(this, p, LocalDateTime.now());
        if (s == null) return;

        long now = SystemClock.elapsedRealtime();
        if (p.equals(lastBlockedPackage) && now - lastBlockedAt < 1500) return;
        lastBlockedPackage = p;
        lastBlockedAt = now;
        foregroundPackage = null;

        performGlobalAction(GLOBAL_ACTION_HOME);
        Intent i = new Intent(this, AsleepActivity.class)
                .addFlags(Intent.FLAG_ACTIVITY_NEW_TASK | Intent.FLAG_ACTIVITY_CLEAR_TASK)
                .putExtra(AsleepActivity.EXTRA_PACKAGE, p)
                .putExtra(AsleepActivity.EXTRA_SCHEDULE_ID, s.id);
        startActivity(i);
    }

    /** Our own app and the home screen are never put to sleep. */
    private void refreshPassThrough() {
        passThrough.clear();
        passThrough.add(getPackageName());
        Intent home = new Intent(Intent.ACTION_MAIN).addCategory(Intent.CATEGORY_HOME);
        List<ResolveInfo> homes = getPackageManager()
                .queryIntentActivities(home, PackageManager.MATCH_DEFAULT_ONLY);
        for (ResolveInfo r : homes) passThrough.add(r.activityInfo.packageName);
    }

    private boolean isInputMethod(String p) {
        InputMethodManager imm = getSystemService(InputMethodManager.class);
        if (imm == null) return false;
        for (InputMethodInfo info : imm.getEnabledInputMethodList()) {
            if (info.getPackageName().equals(p)) return true;
        }
        return false;
    }

    @Override
    public void onInterrupt() {}

    @Override
    public void onDestroy() {
        handler.removeCallbacks(recheck);
        super.onDestroy();
    }
}
